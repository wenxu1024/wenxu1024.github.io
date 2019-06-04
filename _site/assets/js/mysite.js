$(document).ready(function () {

	 $("#sidebar-nav").mCustomScrollbar({
         theme: "minimal"
    });

    $('#sidebarCollapse').on('click', function () {
        $('#sidebar-nav').toggleClass('active');
        $('#content').toggleClass('active');
        $(this).toggleClass('active');
      });
    });
