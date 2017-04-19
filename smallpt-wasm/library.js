mergeInto(LibraryManager.library, {
  getSpp: function() {
  	var name = "spp";
    var reg = new RegExp("(^|&)" + name + "=([^&]*)(&|$)", "i"); 
		var r = window.location.search.substr(1).match(reg); 
		console.log(r);
		if (r != null) return unescape(r[2]); return 16; 
  },
	showPercentage: function(msg) {
		console.log("Rendering: " + msg.toFixed(2) + "%");
	}
});