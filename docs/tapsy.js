// Shared TAPSy widget. Embedded illustrations use their parent page's assistant.
(function(){
 if(window.self !== window.top || document.getElementById('tapsy-widget'))return;
 var root=document.createElement('div');
 root.id='tapsy-widget';
 root.innerHTML="<div style=\"position:fixed;bottom:20px;right:20px;z-index:99999;display:flex;flex-direction:column;align-items:center;gap:4px\">\n  <div id=\"srna-bubble\" style=\"background:#2C5F8A;color:white;font-family:Outfit,sans-serif;font-size:0.78rem;font-weight:600;padding:5px 12px;border-radius:20px;box-shadow:0 2px 8px rgba(44,95,138,0.3);white-space:nowrap;animation:bubble-pulse 4s ease-in-out infinite;\">Need help from TAPSy?</div>\n  \n  <button id=\"srna-chat-btn\" type=\"button\" aria-controls=\"srna-chat-window\" aria-expanded=\"false\" title=\"Ask TAPSy\">\n<img src=\"robot_icon.svg\" alt=\"TAPSy assistant\" style=\"width:80px;height:80px;object-fit:contain;filter:drop-shadow(0 2px 4px rgba(44,95,138,0.2));\">\n</button>\n</div>\n<div id=\"srna-chat-window\" role=\"dialog\" aria-label=\"TAPSy assistant\" inert>\n  <div class=\"chat-header\">\n    <div class=\"chat-header-avatar\" style=\"background:none;padding:2px;\"><img src=\"sRNA-taps_logo_centre.svg\" alt=\"sRNA-TAPS\" style=\"width:32px;height:32px;object-fit:contain;\"></div>\n    <div class=\"chat-header-info\">\n      <div class=\"chat-header-name\">TAPSy</div>\n      <div class=\"chat-header-status\">&#9679; Online</div>\n    </div>\n    <button class=\"chat-close\" type=\"button\" aria-label=\"Close TAPSy\">&#10005;</button>\n  </div>\n  <div class=\"chat-messages\" id=\"srna-chat-messages\" role=\"log\" aria-live=\"polite\">\n    <div class=\"chat-msg bot\"><div class=\"chat-avatar\" style=\"background:none;padding:0;\"><img src=\"robot_icon.svg\" style=\"width:28px;height:28px;object-fit:contain;\"></div><div class=\"chat-bubble\">Hi! I am TAPSy, your sRNA-TAPS assistant. Ask me anything about the pipeline, RNA methylation, or general bioinformatics!</div></div>\n  </div>\n  <div class=\"chat-input-row\">\n    <input type=\"text\" id=\"srna-chat-input\" aria-label=\"Message TAPSy\" placeholder=\"Ask a question...\">\n    <button id=\"srna-chat-send\" type=\"button\" aria-label=\"Send message\"><svg viewBox=\"0 0 24 24\"><path d=\"M2.01 21L23 12 2.01 3 2 10l15 2-15 2z\"/></svg></button>\n  </div>\n</div>";
 document.body.appendChild(root);

var srnaOpen=false;
var srnaHistory=[];
function toggleChat(){
 srnaOpen=!srnaOpen;
 var panel=document.getElementById('srna-chat-window');
 panel.classList.toggle('open',srnaOpen);
 panel.inert=!srnaOpen;
 document.getElementById('srna-chat-btn').setAttribute('aria-expanded',String(srnaOpen));
 document.getElementById('srna-bubble').style.display=srnaOpen?'none':'block';
 (srnaOpen?document.getElementById('srna-chat-input'):document.getElementById('srna-chat-btn')).focus();
}
var pending=false;
async function srnaChat(){
  var input=document.getElementById('srna-chat-input');
  var msg=input.value.trim();
  if(!msg || pending)return;
  pending=true;
  input.value='';
  srnaAddMsg('user',msg);
  srnaHistory.push({role:'user',content:msg});
  var t=srnaTyping();
  try{
    var r=await fetch('https://srna-taps-chat.bennett-henzeler.workers.dev',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({model:'claude-sonnet-4-5',max_tokens:1000,system:'You are TAPSy, the sRNA-TAPS assistant. You are an expert in the sRNA-TAPS pipeline for detecting RNA m5C using TAPS chemistry, small RNA biology (miRNA, tRNA, rRNA, snoRNA), RNA epigenetics, and general bioinformatics. Be concise and helpful.',messages:srnaHistory})});
    var txt=await r.text();
    
    t.remove();
    var d=JSON.parse(txt);
    var reply=(d.content&&d.content[0]&&d.content[0].text)||(d.error&&d.error.message)||('RAW: '+txt);
    srnaHistory.push({role:'assistant',content:reply});
    srnaAddMsg('bot',reply);
  }catch(e){t.remove();srnaAddMsg('bot','Error: '+e.message);}finally{pending=false;}
}
function srnaAddMsg(role,text){
  var m=document.getElementById('srna-chat-messages');
  var d=document.createElement('div');
  d.className='chat-msg '+role;
  var av=role==='bot'?'<div class="chat-avatar" style="background:none;padding:0;"><img src="robot_icon.svg" style="width:28px;height:28px;object-fit:contain;"></div>':'';
  var escaped=String(text).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
  var f=escaped.replace(/`([^`]+)`/g,'<code style="background:#f0eeea;padding:1px 4px;border-radius:3px;font-family:monospace">$1</code>').replace(/\n/g,'<br>');
  d.innerHTML=av+'<div class="chat-bubble">'+f+'</div>';
  m.appendChild(d);
  m.scrollTop=m.scrollHeight;
  return d;
}
function srnaTyping(){
  var m=document.getElementById('srna-chat-messages');
  var d=document.createElement('div');
  d.className='chat-msg bot';
  d.innerHTML='<div class="chat-avatar" style="background:none;padding:0;"><img src="robot_icon.svg" style="width:28px;height:28px;object-fit:contain;"></div><div class="chat-typing"><span></span><span></span><span></span></div>';
  m.appendChild(d);
  m.scrollTop=m.scrollHeight;
  return d;
}

 document.getElementById('srna-chat-btn').addEventListener('click',toggleChat);
 root.querySelector('.chat-close').addEventListener('click',toggleChat);
 document.getElementById('srna-chat-send').addEventListener('click',srnaChat);
 document.getElementById('srna-chat-input').addEventListener('keydown',function(e){if(e.key==='Enter' && !e.isComposing){e.preventDefault();srnaChat();}});
 root.addEventListener('keydown',function(e){if(e.key==='Escape' && srnaOpen){e.stopPropagation();toggleChat();}});
})();
