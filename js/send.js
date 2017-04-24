$(document).ready(function () {
    // Forward unique URL to parent maizegdb.website
    var URL = $(location).attr('href');
    //alert(URL.slice(49));
    var parent = $("#parent-div").html();
    console.log(parent);
    if (parent != "true") {
        window.parent.location.href = "http://www.maizegdb.org/snpversity/send/" + URL.slice(49);
        parent = "false";
    }

//load
    $('#taxaTable').load("html/taxa_codes_modal.html").hide();
    $('#nucleotideTable').load("html/nucleotide_modal.html").hide();

//initialize loading-display for population  of table
    $(document).ajaxStart(displayLoading).ajaxStop($.unblockUI());
//Initialize table:
    rePopulate();

//Slider for zooming in & out:
    $("#slider").slider({
        value: 1,
        min: 0.4,
        max: 1.01,
        step: 0.05,
        slide: function (event, ui) {
            $("#amount").val(ui.value);
            $('#tableOfResults').css({"-moz-transform": "scale(" + ui.value + ")"});
            $('#tableOfResults').css({"-moz-transform-origin": "0% 0%"});
            $('#tableOfResults').css({"zoom": ui.value});
            //$('body,html').animate({scrollTop: 0});
        }
    });
    $("#amount").val($("#slider").slider("value"));  //update value of slider description

    // Convert to CSV option:
    $('.convertCSV').click(
        function() {
            var filename = $("#query-id").html() + ".csv";
            exportTableToCSV.apply(this, [$('#tableOfResults'), filename]);
        });
        
});

function rePopulate() {
    displayLoading();
    var file = $('#pages').find(':selected').val();
    var assembly = $('#version').html();
    var query = $('#query-id').html();
    if (typeof file === 'undefined') {  //For initialisation
        file = $('#first-page').html();
    }

    $.ajax({
        type: "POST",
        async: true,
        url: "get_table_body.php",
        data: 'file=' + JSON.stringify(file) + '&assembly=' + JSON.stringify(assembly), //+'&next_file='+JSON.stringify(next_file),
        beforeSend: function () {
        },
        success: function (data) {
            $('#tableOfResults tbody').html(data);
            if (data == "ERROR"){
                window.open("error.php?query="+query, "_self");
            }
            $('body,html').animate({scrollTop: 0});
            $.unblockUI();
        }
    });
    /* Show table co-ordinates?
    $('#tableOfResults td').click(function () {
        alert('My position in table is: ' + this.cellIndex + 'x' + this.parentNode.rowIndex);
    });*/
}

function togglePageElementVisibility(ele) {
    var obj = typeof ele == 'object'
        ? ele : document.getElementById(ele);

    if (obj.style.display == 'none')
        obj.style.display = 'block';
    else
        obj.style.display = 'none';
    return false;
}

/* Approximate zoom by padding $(function() {
 $("#slider").slider({
 value: 5,
 min: 0,
 max: 5,
 step: 1,
 slide: function (event, ui) {
 $("#amount").val(ui.value);
 $('td').css({"padding": ui.value});
 $('th').css({"padding": ui.value});
 }
 });
 $("#amount").val($("#slider").slider("value"));
 });*/

function showTaxaCodes() {
    $('#taxaTable').draggable();
    $('#taxaTable').click(function () {
    });
    if ($('#showCodeT').text() === "Show Stock Box") {
        $('#showCodeT').text("Hide Stock Box");
        $('#taxaTable').show();
        $('#taxaTable').css({"right": "0", "top": "100px"});
    }
    else {
        $('#showCodeT').text("Show Stock Box");
        $('#taxaTable').hide();
    }
}

function showNucleotideCodes() {
    $('#nucleotideTable').draggable();
    $('#nucleotideTable').click(function () {
        //onclick event
    });
    if ($('#showCodeN').text() === "Show Nucleotide Box") {
        $('#showCodeN').text("Hide Nucleotide Box");
        $('#nucleotideTable').show();
        $('#nucleotideTable').css({"left": "0", "top": "100px"});
    }
    else {
        $('#showCodeN').text("Show Nucleotide Box");
        $('#nucleotideTable').hide();
    }
}
function nextPg() {
    var current = $("#pages").val();
    var countTotal = 0;
    var index = 0;
    $("#pages option").each(function (i) {
        if ($(this).val() == current) {
            index = countTotal;
        }
        else {
            countTotal = countTotal + 1;
        }
    });
    if (index < countTotal) {
        $("#pages option").eq(index + 1).prop("selected", true).change();
    }
}
function prevPg() {
    var current = $("#pages").val();
    var countTotal = 0;
    var index = 0;
    $("#pages option").each(function (i) {
        if ($(this).val() == current) {
            index = countTotal;
        }
        else {
            countTotal = countTotal + 1;
        }
    });
    if (index > 0) {
        $("#pages option").eq(index - 1).prop("selected", true).change();
    }
}

function displayLoading() {
    //checkSubmission();
    $.blockUI({
        message: $("#loading"),
        css: {
            border: 'none',
            padding: '15px',
            backgroundColor: '#000',
            '-webkit-border-radius': '7px',
            '-moz-border-radius': '7px',
            opacity: .5,
            color: '#fff'
        }
    });
}

function exportTableToCSV($table, filename) {
    //var $header = $table.find('th:has(div)');
    var $rows = $.merge($table.find('tr:has(th)'),$table.find('tr:has(td)')) ,

        // Temporary delimiter characters unlikely to be typed by keyboard
        // This is to avoid accidentally splitting the actual contents
        tmpColDelim = String.fromCharCode(11), // vertical tab character
        tmpRowDelim = String.fromCharCode(0), // null character

        // actual delimiter characters for CSV format
        colDelim = '","',
        rowDelim = '"\r\n"',

        // Grab text from table into CSV formatted string
        csv = '"' + $rows.map(function (i, row) {
                var $row = $(row),
                    $cols = $.merge($row.find('div'),$row.find('td'));//$.merge($row.find('div'),($.merge($row.find('td'),$row.find('th'))));

                return $cols.map(function (j, col) {
                    var $col = $(col),
                        text = $col.text();

                    return text.replace('"', '""'); // escape double quotes

                }).get().join(tmpColDelim);

            }).get().join(tmpRowDelim)
                .split(tmpRowDelim).join(rowDelim)
                .split(tmpColDelim).join(colDelim) + '"',

        // Data URI
        csvData = 'data:application/csv;charset=utf-8,' + encodeURIComponent(csv);

    $(this)
        .attr({
            'download': filename,
            'href': csvData,
            'target': '_blank'
        });
}

function moveScroll() {
    var scroll_top = $(window).scrollTop();
    var scroll_left = $(window).scrollLeft();
    var anchor_top = $("#tableOfResults").offset().top;
    var anchor_left = $("#tableOfResults").offset().left;
    var anchor_bottom = $("#bottom_anchor").offset().top;
    
    $("#clone").find("thead").css({
        width: $("#tableOfResults thead").width()+"px",
        position: 'absolute',
        left: - scroll_left  + 'px'
    });
    
    $("#tableOfResults").find(".first").css({
        position: 'absolute',
        left: scroll_left + anchor_left + 'px'
    });
    
    if (scroll_top >= anchor_top && scroll_top <= anchor_bottom) {
        clone_table = $("#clone");
        if (clone_table.length == 0) {
            clone_table = $("#tableOfResults")
                .clone()
                .attr('id', 'clone')
                .css({
                    width: $("#tableOfResults").width()+"px",
                    position: 'fixed',
                    pointerEvents: 'none',
                    left: $("#tableOfResults").offset().left+'px',
                    top: 0
                })
                .appendTo($("#table_container"))
                .css({
                    visibility: 'hidden'
                })
                .find("thead").css({
                    visibility: 'visible'
                });
        }
    }
    else {
        $("#clone").remove();
    }
}
