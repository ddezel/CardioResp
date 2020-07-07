$(document).on('shiny:connected', function(event) {
  var screenInfo = $(window).width();
  if (screenInfo > 550) {
    $("#navbar").toggle();
  }
  //console.log(screenInfo);
  
  $(window).on('resize', function(){
      var win = $(this); //this = window
      if (win.width() > 550) {
        $("#navbar").toggle();
      }
  });
});