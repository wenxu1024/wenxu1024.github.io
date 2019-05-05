$(document).ready(function () {
    $('#sidebarCollapse').on('click', function () {
        $('#sidebar-wrapper').toggleClass('active');
        $(this).toggleClass('active');
      });
    });
