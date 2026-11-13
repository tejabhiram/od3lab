/* OD3 host gene expression atlas. Data precomputed from Maria Taher Plumber's gene_results
   with the same statistics as her tool.ipynb (DESeq2 per study, DerSimonian-Laird pooling,
   Bonferroni across pathogens). */
(function(){
  var BASE = "data/atlas/";
  var PALETTE = ["#4C7FA8","#E3A35C","#6BA368","#B5668F","#8C6BB1"], GREY = "#B8BEC6", DOT = "#6B7280";
  var META = null, GENES = null, UPPER = null, CACHE = {};
  var $ = function(id){ return document.getElementById(id); };
  var NS = "http://www.w3.org/2000/svg";

  function esc(s){ return String(s).replace(/[&<>"]/g, function(c){ return {"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]; }); }
  function disp(p){ return (META.display[p]) || p; }
  function stars(p){ if(p==null) return "ns"; if(p<1e-4) return "****"; if(p<1e-3) return "***"; if(p<1e-2) return "**"; if(p<5e-2) return "*"; return "ns"; }
  function fmt(x,d){ return x==null || !isFinite(x) ? "" : (+x).toFixed(d==null?2:d); }
  function fmtp(x){ if(x==null) return ""; if(x===0) return "< 1e-300"; return x < 0.001 ? (+x).toExponential(1) : (+x).toPrecision(2); }
  function quant(a,q){ var s=a.slice().sort(function(x,y){return x-y}); var pos=(s.length-1)*q, lo=Math.floor(pos), hi=Math.ceil(pos); return s[lo]+(s[hi]-s[lo])*(pos-lo); }
  function mean(a){ return a.reduce(function(s,v){return s+v},0)/a.length; }
  function sd(a){ var m=mean(a); return Math.sqrt(a.reduce(function(s,v){return s+(v-m)*(v-m)},0)/(a.length-1)); }

  // Gaussian KDE with Scott's bandwidth, evaluated over [min, max] like matplotlib's violinplot
  function kde(y){
    var n=y.length, s=sd(y); if(!(s>0)) return null;
    var bw=s*Math.pow(n,-0.2), lo=Math.min.apply(null,y), hi=Math.max.apply(null,y), pts=[], i, j;
    for(i=0;i<100;i++){ var x=lo+(hi-lo)*i/99, d=0; for(j=0;j<n;j++){ var z=(x-y[j])/bw; d+=Math.exp(-0.5*z*z); } pts.push([x,d/(n*bw*Math.sqrt(2*Math.PI))]); }
    return pts;
  }
  function seeded(seed){ return function(){ seed=(seed*16807)%2147483647; return (seed-1)/2147483646; }; }

  function el(tag, attrs, parent){ var e=document.createElementNS(NS,tag); for(var k in attrs) e.setAttribute(k,attrs[k]); if(parent) parent.appendChild(e); return e; }
  function txt(parent, x, y, s, attrs){ var t=el("text", Object.assign({x:x,y:y},attrs||{}), parent); t.textContent=s; return t; }

  function niceTicks(lo, hi, n){
    var span=hi-lo||1, step=Math.pow(10,Math.floor(Math.log10(span/n))), err=span/n/step;
    step*= err>=7.5?10: err>=3.5?5: err>=1.5?2:1;
    var out=[], t=Math.ceil(lo/step)*step; for(;t<=hi+1e-9;t+=step) out.push(+t.toFixed(6)); return out;
  }

  function loadJSON(path){ return fetch(path).then(function(r){ if(!r.ok) throw new Error(r.status); return r.json(); }); }
  function loadGene(g){ if(!CACHE[g]) CACHE[g]=loadJSON(BASE+"g/"+encodeURIComponent(g)+".json"); return CACHE[g]; }

  // ---------- gene search ----------
  function resolve(text){ return UPPER[String(text||"").trim().toUpperCase()] || null; }
  function suggestions(q, k){
    q=q.trim().toUpperCase(); if(!q) return [];
    var pre=[], sub=[];
    for(var i=0;i<GENES.length && pre.length<k;i++){ var u=GENES[i].toUpperCase(); if(u.indexOf(q)===0) pre.push(GENES[i]); }
    if(pre.length<k) for(i=0;i<GENES.length && pre.length+sub.length<k;i++){ u=GENES[i].toUpperCase(); if(u.indexOf(q)>0) sub.push(GENES[i]); }
    return pre.concat(sub);
  }
  function close(q){ // "did you mean" by edit distance on same-length-ish names
    q=q.trim().toUpperCase(); var best=[];
    function lev(a,b){ var d=[],i,j; for(i=0;i<=a.length;i++){d[i]=[i];} for(j=1;j<=b.length;j++)d[0][j]=j;
      for(i=1;i<=a.length;i++)for(j=1;j<=b.length;j++)d[i][j]=Math.min(d[i-1][j]+1,d[i][j-1]+1,d[i-1][j-1]+(a[i-1]===b[j-1]?0:1)); return d[a.length][b.length]; }
    for(var i=0;i<GENES.length;i++){ var u=GENES[i].toUpperCase(); if(Math.abs(u.length-q.length)>2) continue; var d=lev(q,u); if(d<=2) best.push([d,GENES[i]]); }
    return best.sort(function(a,b){return a[0]-b[0]}).slice(0,5).map(function(x){return x[1]});
  }
  function attachAutocomplete(input, list, onPick, multi){
    var active=-1;
    function current(){ if(!multi) return input.value; var parts=input.value.split(/[,\s]+/); return parts[parts.length-1]; }
    function render(){
      var items=suggestions(current(), 8); active=-1;
      list.innerHTML=items.map(function(g){ return '<li role="option" data-g="'+esc(g)+'">'+esc(g)+'</li>'; }).join("");
      list.hidden=!items.length;
    }
    function pick(g){
      if(multi){ var parts=input.value.split(/([,\s]+)/); parts[parts.length-1]=g; input.value=parts.join("")+", "; }
      else input.value=g;
      list.hidden=true; if(!multi) onPick(g); input.focus();
    }
    input.addEventListener("input", render);
    input.addEventListener("keydown", function(e){
      var lis=list.querySelectorAll("li");
      if(e.key==="ArrowDown"&&lis.length){ active=Math.min(active+1,lis.length-1); e.preventDefault(); }
      else if(e.key==="ArrowUp"&&lis.length){ active=Math.max(active-1,0); e.preventDefault(); }
      else if(e.key==="Enter"){ if(active>=0&&!list.hidden){ e.preventDefault(); pick(lis[active].dataset.g); return; } list.hidden=true; return; }
      else if(e.key==="Escape"){ list.hidden=true; return; } else return;
      [].forEach.call(lis,function(li,i){ li.classList.toggle("on",i===active); if(i===active) li.scrollIntoView({block:"nearest"}); });
    });
    list.addEventListener("mousedown", function(e){ var li=e.target.closest("li"); if(li){ e.preventDefault(); pick(li.dataset.g); } });
    input.addEventListener("blur", function(){ setTimeout(function(){ list.hidden=true; },150); });
  }

  // ---------- main violin plot ----------
  function orderOf(rec){ return META.order.filter(function(p){ return rec.p[p]; }).concat(Object.keys(rec.p).filter(function(p){ return META.order.indexOf(p)<0; }).sort()); }
  function infValues(rec, p){
    var out=[]; rec.s.forEach(function(s){ var st=META.studies[s[0]]; if(st.pathogen!==p) return; st.conds.forEach(function(c,i){ if(c==="i" && s[6][i]!=null) out.push(s[6][i]); }); });
    return out;
  }
  function drawViolin(rec, svg, onPathogen){
    svg.innerHTML="";
    var order=orderOf(rec), groups={}, all=[];
    order.forEach(function(p){ groups[p]=infValues(rec,p); all=all.concat(groups[p]); });
    var sig=order.filter(function(p){ return rec.p[p][5]!=null && rec.p[p][5]<0.05; }), colour={};
    sig.forEach(function(p,i){ colour[p]=PALETTE[i%PALETTE.length]; });
    var W=760, H=520, m={l:64,r:16,t:74,b:118}, pw=W-m.l-m.r, ph=H-m.t-m.b;
    svg.setAttribute("viewBox","0 0 "+W+" "+H);
    var lo=Math.min(0,Math.min.apply(null,all)), hi=Math.max(0,Math.max.apply(null,all)), pad=(hi-lo)*0.05||1; lo-=pad; hi+=pad;
    var X=function(i){ return m.l+pw*(i+0.5)/order.length; }, Y=function(v){ return m.t+ph*(hi-v)/(hi-lo); };
    var bg=el("g",{},svg);
    niceTicks(lo,hi,6).forEach(function(t){
      el("line",{x1:m.l,x2:m.l+pw,y1:Y(t),y2:Y(t),stroke:"#E5E7EB","stroke-width":0.8},bg);
      txt(bg,m.l-8,Y(t)+4,String(t),{"text-anchor":"end","font-size":13,fill:"#111"});
    });
    el("line",{x1:m.l,x2:m.l+pw,y1:Y(0),y2:Y(0),stroke:"grey","stroke-dasharray":"4 3","stroke-width":1},bg);
    el("line",{x1:m.l,x2:m.l,y1:m.t,y2:m.t+ph,stroke:"#111"},bg);
    el("line",{x1:m.l,x2:m.l+pw,y1:m.t+ph,y2:m.t+ph,stroke:"#111"},bg);
    var yl=txt(bg,18,m.t+ph/2,"",{"text-anchor":"middle","font-size":14,transform:"rotate(-90 18 "+(m.t+ph/2)+")"});
    yl.innerHTML='log<tspan baseline-shift="sub" font-size="10">2</tspan>FC vs mock';
    txt(svg,W/2,26,rec.g,{"text-anchor":"middle","font-size":19,"font-weight":600,fill:"#0c1b2a"});
    var rng=seeded(7), halfw=0.4*pw/order.length;
    order.forEach(function(p,i){
      var y=groups[p], col=colour[p]||GREY, g=el("g",{class:"vgrp",tabindex:0,role:"button","aria-label":disp(p)+" details"},svg), x=X(i);
      el("rect",{x:x-pw/order.length/2,y:m.t,width:pw/order.length,height:ph,fill:"transparent"},g);
      if(y.length>1){ var d=kde(y); if(d){ var mx=Math.max.apply(null,d.map(function(q){return q[1]}));
        var path=d.map(function(q,k){ return (k?"L":"M")+(x+halfw*q[1]/mx).toFixed(1)+","+Y(q[0]).toFixed(1); }).join("")+
          d.slice().reverse().map(function(q){ return "L"+(x-halfw*q[1]/mx).toFixed(1)+","+Y(q[0]).toFixed(1); }).join("")+"Z";
        el("path",{d:path,fill:col,"fill-opacity":0.45,stroke:col,"stroke-width":1},g); } }
      y.forEach(function(v){ el("circle",{cx:x+(rng()*0.24-0.12)*pw/order.length*1.25,cy:Y(v),r:2.4,fill:colour[p]||DOT,"fill-opacity":0.85},g); });
      if(y.length){ var q1=quant(y,0.25), md=quant(y,0.5), q3=quant(y,0.75), u=pw/order.length;
        el("line",{x1:x-0.2*u,x2:x+0.2*u,y1:Y(md),y2:Y(md),stroke:"#000","stroke-width":2.6},g);
        [q1,q3].forEach(function(q){ el("line",{x1:x-0.1*u,x2:x+0.1*u,y1:Y(q),y2:Y(q),stroke:"#000","stroke-width":1,"stroke-dasharray":"3 2"},g); }); }
      var lab=stars(rec.p[p][5]);
      txt(svg,x,m.t-14,lab,{"text-anchor":"middle","font-size":13,fill:lab==="ns"?"#888":"#000"});
      txt(g,x,m.t+ph+12,disp(p),{"text-anchor":"end","font-size":13,transform:"rotate(-90 "+x+" "+(m.t+ph+12)+")",dy:"0.35em",class:"xlab"});
      var title=el("title",{},g); title.textContent=disp(p)+": n="+rec.p[p][0]+" samples, "+rec.p[p][1]+" studies, pooled log2FC "+fmt(rec.p[p][3])+", padj "+fmtp(rec.p[p][5])+". Click for study details.";
      g.addEventListener("click",function(){ onPathogen(p); });
      g.addEventListener("keydown",function(e){ if(e.key==="Enter"||e.key===" "){ e.preventDefault(); onPathogen(p); } });
    });
  }

  function summaryTable(rec){
    var rows=orderOf(rec).map(function(p){ return [p].concat(rec.p[p]); }).sort(function(a,b){ return (b[4]==null?-1e9:b[4])-(a[4]==null?-1e9:a[4]); });
    return '<table class="atable"><thead><tr><th>Pathogen</th><th>Samples</th><th>Studies</th><th>Median log2FC</th><th>Pooled log2FC</th><th>SE</th><th>Adj. p vs mock</th><th></th></tr></thead><tbody>'+
      rows.map(function(r){ return '<tr data-p="'+esc(r[0])+'"><td><button class="linkbtn">'+esc(disp(r[0]))+'</button></td><td>'+r[1]+'</td><td>'+r[2]+'</td><td>'+fmt(r[3])+'</td><td><b>'+fmt(r[4])+'</b></td><td>'+fmt(r[5])+'</td><td>'+fmtp(r[6])+'</td><td>'+stars(r[6])+'</td></tr>'; }).join("")+'</tbody></table>';
  }

  // ---------- per-pathogen detail ----------
  function drawDetail(rec, p, box){
    var rows=rec.s.filter(function(s){ return META.studies[s[0]].pathogen===p; });
    var W=760, H=330, m={l:58,r:14,t:34,b:96};
    function panel(title, ylab, build){
      var svg=document.createElementNS(NS,"svg"); svg.setAttribute("viewBox","0 0 "+W+" "+H); svg.setAttribute("class","dplot");
      txt(svg,W/2,18,title,{"text-anchor":"middle","font-size":14,"font-weight":600});
      var yl=txt(svg,14,m.t+(H-m.t-m.b)/2,ylab,{"text-anchor":"middle","font-size":12,transform:"rotate(-90 14 "+(m.t+(H-m.t-m.b)/2)+")"});
      build(svg); return svg;
    }
    var pw=W-m.l-m.r, ph=H-m.t-m.b, n=rows.length, X=function(i){ return m.l+pw*(i+0.5)/n; };
    function axes(svg, lo, hi, Y, zero){
      niceTicks(lo,hi,5).forEach(function(t){ el("line",{x1:m.l,x2:m.l+pw,y1:Y(t),y2:Y(t),stroke:"#EEF0F2"},svg); txt(svg,m.l-6,Y(t)+4,String(t),{"text-anchor":"end","font-size":11}); });
      if(zero) el("line",{x1:m.l,x2:m.l+pw,y1:Y(0),y2:Y(0),stroke:"grey","stroke-dasharray":"4 3"},svg);
      el("line",{x1:m.l,x2:m.l,y1:m.t,y2:m.t+ph,stroke:"#333"},svg); el("line",{x1:m.l,x2:m.l+pw,y1:m.t+ph,y2:m.t+ph,stroke:"#333"},svg);
      rows.forEach(function(s,i){ var st=META.studies[s[0]].study; txt(svg,X(i),m.t+ph+8,st,{"text-anchor":"end","font-size":10.5,transform:"rotate(-45 "+X(i)+" "+(m.t+ph+8)+")",dy:"0.35em"}); });
    }
    // TPM by condition
    var tp=[]; rows.forEach(function(s){ s[5].forEach(function(v){ if(v!=null) tp.push(v); }); });
    var tHi=(Math.max.apply(null,tp.concat([1])))*1.08, Yt=function(v){ return m.t+ph*(tHi-v)/tHi; };
    var rng=seeded(11);
    var a=panel(disp(p)+": "+rec.g+" expression by study","TPM",function(svg){
      axes(svg,0,tHi,Yt,false);
      rows.forEach(function(s,i){ var st=META.studies[s[0]], u=pw/n;
        ["m","i"].forEach(function(c,k){ var vals=[]; st.conds.forEach(function(cc,j){ if(cc===c && s[5][j]!=null) vals.push(s[5][j]); });
          var x=X(i)+(k?0.18:-0.18)*u, fill=k?"#d6604d":"#9e9e9e", dot=k?"#67001f":"#424242";
          if(vals.length>1){ var q1=quant(vals,.25),q3=quant(vals,.75),md=quant(vals,.5),bw=Math.min(0.3*u,22);
            el("rect",{x:x-bw/2,y:Yt(q3),width:bw,height:Math.max(1,Yt(q1)-Yt(q3)),fill:fill,"fill-opacity":0.55,stroke:"#555","stroke-width":0.8},svg);
            el("line",{x1:x-bw/2,x2:x+bw/2,y1:Yt(md),y2:Yt(md),stroke:"#222","stroke-width":1.4},svg); }
          vals.forEach(function(v){ el("circle",{cx:x+(rng()-0.5)*Math.min(0.2*u,14),cy:Yt(v),r:2.6,fill:dot},svg); });
        });
      });
      el("rect",{x:W-170,y:m.t-26,width:10,height:10,fill:"#9e9e9e"},svg); txt(svg,W-155,m.t-17,"mock",{"font-size":11});
      el("rect",{x:W-110,y:m.t-26,width:10,height:10,fill:"#d6604d"},svg); txt(svg,W-95,m.t-17,"infected",{"font-size":11});
    });
    // log2FC per study with DESeq2 estimate
    var lv=[]; rows.forEach(function(s){ var st=META.studies[s[0]]; st.conds.forEach(function(c,j){ if(c==="i"&&s[6][j]!=null) lv.push(s[6][j]); });
      if(s[1]!=null){ lv.push(s[1]); if(s[2]!=null){ lv.push(s[1]+1.96*s[2]); lv.push(s[1]-1.96*s[2]); } } });
    var lo=Math.min(0,Math.min.apply(null,lv)), hi=Math.max(0,Math.max.apply(null,lv)), pd=(hi-lo)*0.08||1; lo-=pd; hi+=pd;
    var Yl=function(v){ return m.t+ph*(hi-v)/(hi-lo); };
    var b=panel(disp(p)+": "+rec.g+" log2FC infected vs mock","log2 fold change",function(svg){
      axes(svg,lo,hi,Yl,true);
      rows.forEach(function(s,i){ var st=META.studies[s[0]], x=X(i);
        st.conds.forEach(function(c,j){ if(c==="i"&&s[6][j]!=null) el("circle",{cx:x+(rng()-0.5)*Math.min(0.3*pw/n,18),cy:Yl(s[6][j]),r:3,fill:"#d6604d","fill-opacity":0.85},svg); });
        if(s[1]!=null){
          if(s[2]!=null){ el("line",{x1:x,x2:x,y1:Yl(s[1]+1.96*s[2]),y2:Yl(s[1]-1.96*s[2]),stroke:"#000","stroke-width":1.4},svg);
            [1,-1].forEach(function(sg){ el("line",{x1:x-5,x2:x+5,y1:Yl(s[1]+sg*1.96*s[2]),y2:Yl(s[1]+sg*1.96*s[2]),stroke:"#000","stroke-width":1.4},svg); }); }
          el("path",{d:"M"+x+","+(Yl(s[1])-6)+"L"+(x+6)+","+Yl(s[1])+"L"+x+","+(Yl(s[1])+6)+"L"+(x-6)+","+Yl(s[1])+"Z",fill:s[2]!=null?"#000":"#fff",stroke:"#000"},svg);
        }
        txt(svg,x,m.t-4,META.calls[s[4]],{"text-anchor":"middle","font-size":10.5,fill:"#333"});
      });
    });
    var table='<table class="atable small"><thead><tr><th>Study</th><th>Mock</th><th>Infected</th><th>log2FC</th><th>SE</th><th>padj</th><th>Call</th></tr></thead><tbody>'+
      rows.map(function(s){ var st=META.studies[s[0]]; return '<tr><td><a href="https://www.ncbi.nlm.nih.gov/bioproject/'+esc(st.study)+'" target="_blank" rel="noopener">'+esc(st.study)+'</a></td><td>'+st.n_mock+'</td><td>'+st.n_inf+'</td><td>'+fmt(s[1])+'</td><td>'+fmt(s[2])+'</td><td>'+fmtp(s[3])+'</td><td>'+META.calls[s[4]]+'</td></tr>'; }).join("")+'</tbody></table>';
    box.innerHTML='<div class="detail-head"><h3>'+esc(rec.g)+' in '+esc(disp(p))+'</h3><button class="chip" id="detail-close">Close</button></div>'+
      '<p class="muted small-note">Left: TPM for mock and infected samples in each study. Right: dots are infected samples relative to the study’s mock mean, diamonds are the DESeq2 estimate with 95% CI (open diamond means no replicates). Labels show each study’s call at adjusted p &lt; 0.05.</p>';
    var grid=document.createElement("div"); grid.className="detail-grid"; grid.appendChild(a); grid.appendChild(b); box.appendChild(grid);
    box.insertAdjacentHTML("beforeend",table); box.hidden=false;
    $("detail-close").onclick=function(){ box.hidden=true; };
    box.scrollIntoView({behavior:"smooth",block:"start"});
  }

  // ---------- downloads ----------
  function download(name, text, type){ var a=document.createElement("a"); a.href=URL.createObjectURL(new Blob([text],{type:type||"text/csv"})); a.download=name; document.body.appendChild(a); a.click(); setTimeout(function(){ URL.revokeObjectURL(a.href); a.remove(); },500); }
  function csvRow(a){ return a.map(function(v){ v=v==null?"":String(v); return /[",\n]/.test(v)?'"'+v.replace(/"/g,'""')+'"':v; }).join(","); }
  function samplesCSV(rec){
    var out=[csvRow(["gene","pathogen","study","sample","condition","tpm","log2FC_vs_mock"])];
    rec.s.forEach(function(s){ var st=META.studies[s[0]]; st.samples.forEach(function(id,j){ out.push(csvRow([rec.g,disp(st.pathogen),st.study,id,st.conds[j]==="m"?"mock":"infected",s[5][j],s[6][j]])); }); });
    return out.join("\n");
  }
  function studiesCSV(rec){
    var out=[csvRow(["gene","pathogen","study","n_mock","n_infected","log2FC","lfcSE","padj","call"])];
    rec.s.forEach(function(s){ var st=META.studies[s[0]]; out.push(csvRow([rec.g,disp(st.pathogen),st.study,st.n_mock,st.n_inf,s[1],s[2],s[3],META.calls[s[4]]])); });
    out.push(""); out.push(csvRow(["gene","pathogen","n_samples","n_studies","median_log2FC","pooled_log2FC","pooled_SE","padj_vs_mock"]));
    orderOf(rec).forEach(function(p){ out.push(csvRow([rec.g,disp(p)].concat(rec.p[p]))); });
    return out.join("\n");
  }
  function svgText(svg){ var c=svg.cloneNode(true); c.setAttribute("xmlns",NS); c.setAttribute("style","font-family:Inter,Arial,sans-serif;background:#fff"); return '<?xml version="1.0" encoding="UTF-8"?>\n'+c.outerHTML; }
  function downloadPNG(svg, name){
    var vb=svg.viewBox.baseVal, img=new Image(), s=3;
    img.onload=function(){ var cv=document.createElement("canvas"); cv.width=vb.width*s; cv.height=vb.height*s; var cx=cv.getContext("2d"); cx.fillStyle="#fff"; cx.fillRect(0,0,cv.width,cv.height); cx.drawImage(img,0,0,cv.width,cv.height);
      cv.toBlob(function(b){ var a=document.createElement("a"); a.href=URL.createObjectURL(b); a.download=name; a.click(); }); };
    img.src="data:image/svg+xml;charset=utf-8,"+encodeURIComponent(svgText(svg));
  }

  // ---------- multi-gene heatmap ----------
  function drawHeatmap(recs, svg){
    svg.innerHTML="";
    var paths=META.order.filter(function(p){ return recs.some(function(r){ return r.p[p]; }); });
    var cw=58, ch=28, m={l:110,t:110,r:20,b:50}, W=m.l+cw*paths.length+m.r, H=m.t+ch*recs.length+m.b;
    svg.setAttribute("viewBox","0 0 "+W+" "+H);
    var vmax=0; recs.forEach(function(r){ paths.forEach(function(p){ if(r.p[p]&&r.p[p][3]!=null) vmax=Math.max(vmax,Math.abs(r.p[p][3])); }); }); vmax=Math.max(1,Math.ceil(vmax));
    function colour(v){ if(v==null) return "#f3f4f6"; var t=Math.max(-1,Math.min(1,v/vmax)), c1=t>0?[178,24,43]:[33,102,172], k=Math.abs(t);
      return "rgb("+[0,1,2].map(function(i){ return Math.round(255+(c1[i]-255)*k); }).join(",")+")"; }
    paths.forEach(function(p,j){ var x=m.l+cw*j+cw/2; txt(svg,x,m.t-8,disp(p),{"text-anchor":"start","font-size":12,transform:"rotate(-50 "+x+" "+(m.t-8)+")"}); });
    recs.forEach(function(r,i){ var y=m.t+ch*i;
      txt(svg,m.l-8,y+ch/2+4,r.g,{"text-anchor":"end","font-size":12.5,"font-weight":600});
      paths.forEach(function(p,j){ var v=r.p[p]?r.p[p][3]:null, pa=r.p[p]?r.p[p][5]:null, x=m.l+cw*j;
        var rect=el("rect",{x:x+1,y:y+1,width:cw-2,height:ch-2,fill:colour(v),rx:3},svg);
        var t=el("title",{},rect); t.textContent=r.g+" in "+disp(p)+": pooled log2FC "+fmt(v)+", padj "+fmtp(pa);
        var st=stars(pa); if(st!=="ns") txt(svg,x+cw/2,y+ch/2+4,st,{"text-anchor":"middle","font-size":11,fill:Math.abs(v||0)/vmax>0.55?"#fff":"#111"});
      });
    });
    var ly=m.t+ch*recs.length+18, lx=m.l, lw=Math.min(240,cw*paths.length);
    for(var k=0;k<=40;k++){ var v=-vmax+2*vmax*k/40; el("rect",{x:lx+lw*k/41,y:ly,width:lw/41+0.5,height:10,fill:colour(v)},svg); }
    txt(svg,lx,ly+24,"-"+vmax,{"font-size":11}); txt(svg,lx+lw/2,ly+24,"0",{"font-size":11,"text-anchor":"middle"}); txt(svg,lx+lw,ly+24,"+"+vmax,{"font-size":11,"text-anchor":"end"});
    txt(svg,lx+lw+12,ly+9,"pooled log2FC vs mock (stars: Bonferroni-adjusted p)",{"font-size":11,fill:"#555"});
  }

  // ---------- wiring ----------
  var current=null;
  function showGene(g){
    var hit=resolve(g), msg=$("gene-msg");
    if(!hit){ var c=close(g); msg.innerHTML="<b>"+esc(g)+"</b> is not in the gene pool."+(c.length?" Did you mean "+c.map(function(x){ return '<button class="linkbtn" data-g="'+esc(x)+'">'+esc(x)+'</button>'; }).join(", ")+"?":""); msg.hidden=false;
      [].forEach.call(msg.querySelectorAll("button"),function(b){ b.onclick=function(){ $("gene-input").value=b.dataset.g; showGene(b.dataset.g); }; }); return; }
    msg.hidden=true; $("gene-input").value=hit; $("result").setAttribute("aria-busy","true");
    loadGene(hit).then(function(rec){
      current=rec; $("result").hidden=false; $("result").setAttribute("aria-busy","false");
      drawViolin(rec,$("violin"),function(p){ drawDetail(rec,p,$("detail")); });
      $("summary").innerHTML=summaryTable(rec);
      [].forEach.call($("summary").querySelectorAll("tr[data-p] button"),function(b){ b.onclick=function(){ drawDetail(rec,b.closest("tr").dataset.p,$("detail")); }; });
      $("detail").hidden=true;
      var n=rec.s.length; $("result-note").textContent=rec.g+" was quantified in "+n+" of "+META.studies.length+" studies. Click a pathogen in the plot or table for study-level detail.";
      try{ history.replaceState(null,"","#"+encodeURIComponent(hit)); }catch(e){}
    }).catch(function(){ msg.textContent="Could not load data for "+hit+". Please try again."; msg.hidden=false; });
  }
  function compare(){
    var names=$("multi-input").value.split(/[,\s]+/).filter(Boolean), ok=[], bad=[];
    names.forEach(function(n){ var h=resolve(n); if(h){ if(ok.indexOf(h)<0) ok.push(h); } else bad.push(n); });
    ok=ok.slice(0,30);
    $("multi-msg").textContent=bad.length?"Not found: "+bad.join(", "):""; $("multi-msg").hidden=!bad.length;
    if(!ok.length) return;
    Promise.all(ok.map(loadGene)).then(function(recs){ $("multi-result").hidden=false; drawHeatmap(recs,$("heatmap"));
      $("multi-csv").onclick=function(){ var paths=META.order, out=[csvRow(["gene"].concat(paths.map(function(p){return disp(p)+"_pooled_log2FC"})).concat(paths.map(function(p){return disp(p)+"_padj"})))];
        recs.forEach(function(r){ out.push(csvRow([r.g].concat(paths.map(function(p){return r.p[p]?r.p[p][3]:""})).concat(paths.map(function(p){return r.p[p]?r.p[p][5]:""})))); });
        download("od3_atlas_"+ok.length+"_genes.csv",out.join("\n")); };
      $("multi-png").onclick=function(){ downloadPNG($("heatmap"),"od3_atlas_heatmap.png"); };
    });
  }

  Promise.all([loadJSON(BASE+"meta.json"),loadJSON(BASE+"genes.json")]).then(function(res){
    META=res[0]; GENES=res[1]; UPPER={}; GENES.forEach(function(g){ UPPER[g.toUpperCase()]=g; });
    $("gene-count").textContent=GENES.length.toLocaleString();
    $("study-count").textContent=META.studies.length;
    $("pathogen-count").textContent=META.order.length;
    $("loading").hidden=true; $("tool").hidden=false;
    attachAutocomplete($("gene-input"),$("gene-list"),showGene,false);
    attachAutocomplete($("multi-input"),$("multi-list"),null,true);
    $("gene-form").onsubmit=function(e){ e.preventDefault(); showGene($("gene-input").value); };
    $("multi-form").onsubmit=function(e){ e.preventDefault(); compare(); };
    [].forEach.call(document.querySelectorAll("[data-example]"),function(b){ b.onclick=function(){ showGene(b.dataset.example); }; });
    $("dl-samples").onclick=function(){ if(current) download(current.g+"_all_samples.csv",samplesCSV(current)); };
    $("dl-studies").onclick=function(){ if(current) download(current.g+"_studies_and_pooled.csv",studiesCSV(current)); };
    $("dl-png").onclick=function(){ if(current) downloadPNG($("violin"),current.g+"_violin.png"); };
    $("dl-svg").onclick=function(){ if(current) download(current.g+"_violin.svg",svgText($("violin")),"image/svg+xml"); };
    var h=decodeURIComponent((location.hash||"").slice(1)); if(h && resolve(h)) showGene(h);
  }).catch(function(){ $("loading").textContent="The atlas data could not be loaded. Please refresh the page."; });
})();
