(function(){
  var LOGO='<svg viewBox="0 0 40 40" aria-hidden="true"><circle cx="20" cy="20" r="18.5" fill="none" stroke="#0c1b2a" stroke-width="1.6"/><circle cx="20" cy="20" r="11" fill="#0d6e78"/><circle cx="20" cy="20" r="4.6" fill="#0c1b2a"/><circle cx="23.6" cy="16.4" r="2" fill="#fff"/><path d="M20 1.5v5M20 33.5v5M1.5 20h5M33.5 20h5" stroke="#d8892b" stroke-width="1.6" stroke-linecap="round"/></svg>';
  var pages=[["index.html","Home"],["research.html","Research"],["people.html","People"],["publications.html","Publications"],["news.html","News"],["funding.html","Funding"],["tools.html","Tools"]];
  var here=(location.pathname.split("/").pop()||"index.html");
  var nav=pages.map(function(p){return '<a href="'+p[0]+'"'+(p[0]===here?' class="active" aria-current="page"':'')+'>'+p[1]+'</a>'}).join("")+'<a href="join.html" class="cta'+(here==="join.html"?" active":"")+'">Join / Contact</a>';
  var h=document.getElementById("site-header");
  if(h){h.innerHTML='<div class="wrap"><a class="brand" href="index.html">'+LOGO+'<div><b>OD3 Lab</b><span>Ocular Drug Delivery &amp; Discovery</span></div></a><button class="menu-btn" aria-label="Menu" aria-expanded="false">Menu</button><nav class="nav">'+nav+'</nav></div>';
    var b=h.querySelector(".menu-btn"),n=h.querySelector(".nav");b.addEventListener("click",function(){var o=n.classList.toggle("open");b.setAttribute("aria-expanded",o)});}
  var f=document.getElementById("site-footer");
  if(f){f.innerHTML='<div class="wrap"><div class="cols"><div><h4>OD3 Lab</h4><p>Ocular Drug Delivery and Discovery Lab<br>Department of Ophthalmology and Visual Sciences<br>University of Illinois Chicago</p></div><div><h4>Visit</h4><p>909 S. Wolcott Ave, Room COMRB 1080<br>Chicago, IL 60612</p><p><a href="mailto:yteja@uic.edu">yteja@uic.edu</a></p></div><div><h4>Elsewhere</h4><p><a href="https://scholar.google.com/citations?user=M5iDh0UAAAAJ&hl=en">Google Scholar</a><br><a href="https://orcid.org/0000-0001-8584-2685">ORCID</a><br><a href="https://chicago.medicine.uic.edu/ophthalmology-visual-sciences/">UIC Ophthalmology</a><br><a href="join.html">Open positions</a></p></div></div><div class="base"><span>&copy; '+new Date().getFullYear()+' OD3 Lab, University of Illinois Chicago</span><span>Supported by the National Eye Institute</span></div></div>';}
})();

/* Play section videos while they are on screen, pause when scrolled away, resume on return.
   A video the visitor paused by hand stays paused until they press play again. */
(function(){
  var vids=[].slice.call(document.querySelectorAll('video[data-autoplay]'));
  if(!vids.length||!('IntersectionObserver' in window))return;
  if(window.matchMedia&&window.matchMedia('(prefers-reduced-motion: reduce)').matches)return;
  vids.forEach(function(v){
    v.muted=true;v._auto=false;v._userPaused=false;
    v.addEventListener('pause',function(){if(!v._auto&&!v.ended)v._userPaused=true;v._auto=false});
    v.addEventListener('play',function(){v._userPaused=false});
  });
  var io=new IntersectionObserver(function(entries){
    entries.forEach(function(e){
      var v=e.target;
      if(e.isIntersecting&&e.intersectionRatio>=0.6){
        if(!v._userPaused&&v.paused&&!v.ended){var p=v.play();if(p&&p.catch)p.catch(function(){})}
      }else if(!v.paused){v._auto=true;v.pause()}
    });
  },{threshold:[0,0.6]});
  vids.forEach(function(v){io.observe(v)});
})();
