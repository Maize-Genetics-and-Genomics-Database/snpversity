$(document).ready(function () {
    var tooltips = $("[title]").tooltip({
        position: {
            my: "left top",
            at: "right+5 top-5"
        }
    });

    $.noty.defaults.layout = 'topRight';
    $.noty.defaults.theme = 'relax';
    $("#errorMsg").hide();
    //$("#dataSet").hide();
    $("#resultsMax").val(50);
    $("#model-select").focusout(function(){
        getRangeOfGeneModel(document.getElementById('model-select').value, document.getElementById('offset').value);
    });
});


function validateForm(blockScreen) {
    var file = $('#stockFile');
    var dataSet = $("#dataSet").val();
    var startPos = $("#startPosition").val();
    var endPos = $("#endPosition").val();
    if (!isTaxaFile(file.val()) && !valUndefined(file.val())) {
        file.val("");
        return failValidation(' This custom file is not supported. Please select a file with .csv or .stockinfo extension.');
    }
    else if (valUndefined(file.val()) && valUndefined($("#taxa option:selected").eq(0).val())) {
        return failValidation(' Please select a taxon or a file.');
    }
    else if (dataSet == "HapMapV3" && isCsvFile(file.val())) {  // .CSV doesn't make sense for HapMapV3 because there are no projects!
        return failValidation(".CSV files are not supported for HapMapV3. Please submit a .stockinfo file instead.");
    }
    else if ((parseInt(startPos, 10) > parseInt(endPos, 10))) {
        return failValidation("Please check your start and end positions.")
    }
    else {
        //Success:
        if (blockScreen){
            displayLoading();
        }
        $("#errorMsg").hide();

        return true;
    }
}
function onAssemblyChange() {
    var assembly = $("#assembly").val();
    if (assembly == "v2") {
        var v2html = '<option value="ZeaGBSv27publicImputed20150114">ZeaGBSv27publicImputed20150114</option><option value="AllZeaGBSv27public20140528">AllZeaGBSv27public20140528</option>';
        $("#dataSet").html(v2html);
        $("#taxa").val('').trigger("chosen:updated");  //remove any taxa from old dataset TODO: only if assembly changed.
    }
    else {
        var v3html = '<option value="ZeaHM321_raw">HapMapV3</option><option value="ZeaHM321_LinkImpute">HapMapV3 Imputed</option>';
        $("#dataSet").html(v3html);
    }
    onDataSetChange();
    $("#dataSet").show();
    onChromosomeChange();
    $("#model-select").val('');
}
//Only needs to be run if assembly# changed.
function onDataSetChange() {
    var assembly = $("#assembly").val();
    var dataSet = $("#dataSet").val();
    if (assembly == "v3") {
        //hide project-selection
        $("#projects").hide();
        //remove chromosome 0 in options
        $("#chromosome option:contains('0')").eq(0).prop("disabled", true);
        $("#chromosome").val(1);
        onChromosomeChange();
        //populate taxa:
        $.ajax({
            type: "POST",
            data: {action: "getHapMapLines"},
            url: "get_taxa_hapmapv3.php",
            success: function (data) {
                $("#taxa").html(data);
                $("#taxa").trigger("chosen:updated");
            }
        });
        //show change in co-ord assembly
    }
    else {
        $("#projects").show();
        $("#chromosome option:contains('0')").eq(0).prop("disabled", false);
        //show change in co-ord assembly
    }
}

function estimate() {
    var validated = validateForm(false);
    if (validated) {
        var assembly = $("#assembly").val();
        var n = $('.notify-box').noty({
            text: 'Estimating Time...'
        });
        n.setType('alert');
        // get formData, attach file
        var formData = new FormData($('form')[0]);
        var stockFile = $("#stockFile")[0].files[0];
        formData.append('stockFile', stockFile);
        // Submit to PHP
        $.ajax({
            type: "POST",
            data: formData,
            cache: false,
            processData: false,
            contentType: false,
            url: "get_time_estimate.php",
            success: function (data) {
                n.setType('success');
                var time_obj = JSON.parse(data);
                time_obj.timeInMinutes = Math.ceil(time_obj.time / 60);
                if (time_obj.timeInMinutes > 4 && time_obj.timeInMinutes < 10){
                    n.setType('warning');
                }
                else if (time_obj.timeInMinutes > 9){
                    n.setType('error');
                }
                msg = "Time estimate of query using " + time_obj.stocks + " stock(s) and " + time_obj.range + " bp's is: <strong>" + time_obj.timeInMinutes + " minute(s).</strong>";
                n.setText(msg);
            }
        });
        }
}

function onOutputFormatChange() {
    var out = $("#outputFormat").val();
    if (out == 'json') {
        $("#resultsMax").show();
        $("#resultsMaxLabel").show();
    }
    else {
        $("#resultsMax").hide();
        $("#resultsMaxLabel").hide();
    }
}
function onChromosomeChange() {
    var chr = $("#chromosome").val();
    var min = 0;
    var max = 1;
    var min_end = 1;
    if ($("#positions").val() == 'range' && $("#assembly").val() != "v3") {
        switch (chr) {
            case '0':
                min = 1;
                min_end = 5518;
                max = 6932521;
                break;
            case '1':
                min = 1;
                min_end = 6370;
                max = 301331039;
                break;
            case '2':
                min = 1;
                min_end = 9700;
                max = 237042811;
                break;
            case '3':
                min = 1;
                min_end = 32942;
                max = 232096209;
                break;
            case '4':
                min = 1;
                min_end = 28990;
                max = 241426350;
                break;
            case '5':
                min = 1;
                min_end = 281;
                max = 217748593;
                break;
            case '6':
                min = 1;
                min_end = 119364;
                max = 169132442;
                break;
            case '7':
                min = 1;
                min_end = 27;
                max = 176759698;
                break;
            case '8':
                min = 1;
                min_end = 17932;
                max = 175735562;
                break;
            case '9':
                min = 1;
                min_end = 66658;
                max = 156591838;
                break;
            case '10':
                min = 1;
                min_end = 2918;
                max = 150146791;
                break;
            default:
                min = 1;
                max = 2;
        }
    }
    else {
        switch (chr) {
            case '1':
                min = 1;
                min_end = 10004;
                max = 301410279;
                break;
            case '2':
                min = 1;
                min_end = 10056;
                max = 237798960;
                break;
            case '3':
                min = 1;
                min_end = 32942;
                max = 232184005;
                break;
            case '4':
                min = 1;
                min_end = 29023;
                max = 241982791;
                break;
            case '5':
                min = 1;
                min_end = 313;
                max = 217804155;
                break;
            case '6':
                min = 1;
                min_end = 120018;
                max = 169339845;
                break;
            case '7':
                min = 1;
                min_end = 65;
                max = 176174508;
                break;
            case '8':
                min = 1;
                min_end = 17958;
                max = 175289467;
                break;
            case '9':
                min = 1;
                min_end = 66658;
                max = 156862189;
                break;
            case '10':
                min = 1;
                min_end = 228730;
                max = 149584883;
                break;
            default:
                min = 1;
                max = 301410279;
        }
    }
    $("#startPosition").val(min);
    $("#endPosition").val(max);
    $("#endPosition").attr("min", min_end);
}

function onPositionsChange() {
    var positions = $("#positions").val();
    if (positions == "all") {
        $("#startPosition").hide();
        $("#startPosition").prop("required", false);
        $("#endPosition").hide();
        $("#endPosition").prop("required", false);
        $("#startPositionLabel").hide();
        $("#endPositionLabel").hide();
    }
    if (positions == "range")
    {
        $("#startPosition").show();
        $("#startPosition").prop("required", true);
        $("#endPosition").show();
        $("#endPosition").prop("required", true);
        $("#startPositionLabel").show();
        $("#endPositionLabel").show();
    }
}
//takes two iterations to update
function selectMediumQuery(proj) {
    $("#assembly").val("v2");
    onAssemblyChange();
    $("#chromosome").val("1");
    $("#positions").val("range");
    onPositionsChange();
    $("#startPosition").val("10045");
    $("#endPosition").val("1004497");
    $("#dataSet").val("ZeaGBSv27publicImputed20150114");
    $("#project option:contains('" + proj + "')").prop("selected", true).trigger("chosen:updated");
    var defer = $.Deferred(),
        filtered = defer.then(function () {
            return getGBSTaxa();
        });
    defer.resolve();
    filtered.done(function () {
        selectFirstTaxon(proj);
    });
    $("#outputFormat").val("json");
}
;
function selectBigQuery(proj) {
    $("#assembly").val("v2");
    onAssemblyChange();
    $("#chromosome").val("10");
    $("#positions").val("all");
    onPositionsChange();
    $("#dataSet").val("ZeaGBSv27publicImputed20150114");
    $("#project option:contains('" + proj + "')").prop("selected", true).trigger("chosen:updated");
    var defer = $.Deferred(),
        filtered = defer.then(function () {
            return getGBSTaxa();
        });
    defer.resolve();
    filtered.done(function () {
        selectFirstTaxon(proj);
    });
    $("#outputFormat").val("json");
}

function selectSmallQuery(proj) {
    $("#assembly").val("v2");
    onAssemblyChange();
    $("#chromosome").val("0");
    $("#positions").val("range");
    onPositionsChange();
    $("#dataSet").val("ZeaGBSv27publicImputed20150114");
    onDataSetChange();
    $("#project option:contains('" + proj + "')").prop("selected", true).trigger("chosen:updated");
    $("#outputFormat").val("json");
    var defer = $.Deferred(),
        filtered = defer.then(function () {
            return getGBSTaxa();
        });
    defer.resolve();
    filtered.done(function () {
        selectMultipleGBSTaxa(proj);
    });
    $("#startPosition").val("5518");
    $("#endPosition").val("35161");
}

function selectMultipleGBSTaxa(proj) {
    setTimeout(function () {
        for (var i = 1; i < 10; i++) {
            //timeout needed because otherwise chosen will not have updated...Potentially chain the 2 functions?
            $("#taxa").children('optgroup[label="' + proj + '"]').filter(function () {
                $(this).children('option').eq(i).prop("selected", true).trigger('chosen:updated');
            });
        }
    }, 500);
    $("#taxa").children('option').eq(0).prop("selected", false).trigger('chosen:updated');
}

function selectFirstTaxon(proj) {
    if (proj == 'NAM') {
        $("#resultsMax").val(30);
    }
    setTimeout(function () {  //timeout needed because otherwise chosen will not have updated...Potentially chain the 2 functions?
        $("#taxa").children('optgroup[label="' + proj + '"]').filter(function () {
            $(this).children('option').eq(0).prop("selected", true).trigger('chosen:updated');
        });
    }, 500);
}


function failValidation(msg) {
    $('#errorMsg').text(msg).show();
    alert(msg);
    return false;
}
function valUndefined(val) {
    return (val === undefined || val === null || val === '');
}

function getExtension(filename) {
    var parts = filename.split('.');
    return parts[parts.length - 1];
}

function isTaxaFile(filename) {
    var ext = getExtension(filename);
    switch (ext.toLowerCase()) {
        case 'csv':
        case 'stockinfo':
        case 'taxainfo':
        case 'txt':
            return true;
    }
    return false;
}

function isCsvFile(filename) {
    var ext = getExtension(filename);
    return (ext.toLowerCase() == 'csv');
}

function displayLoading() {
    $.blockUI({message: $("#loading"),
        css: {
            border: 'none',
            padding: '15px',
            backgroundColor: '#000',
            '-webkit-border-radius': '7px',
            '-moz-border-radius': '7px',
            opacity: .5,
            color: '#fff'}
    });
}
function getGBSTaxa() {
    var projs_selected = {};
    var i = 0;
    $("#project option:selected").each(function () {
        if ($(this).html() == "R&amp;D") {
            projs_selected[i] = "RandD";
        }
        else {
            projs_selected[i] = $(this).html();
        }
        i += 1;
    });
    var projs_selected_json = JSON.stringify(projs_selected);
    $.ajax({
        type: "POST",
        url: "get_taxa_allzeagbs.php",
        data: 'projects=' + projs_selected_json,
        success: function (data) {
            var taxa_selected = {};
            var j = 0;
            //preserve selections by storing in taxa_selected
            $("#taxa option:selected").each(function () {
                taxa_selected[j] = $(this).html();
                //TODO: Ensure 'All R&D' works alert(taxa_selected[j]);
                j += 1;
            });
            //Update taxa lines:
            $("#taxa").html(data);
            //Make sure that selections are kept (reselect)
            $.each(taxa_selected, function (index, id_selected) {
                $("#taxa option:contains('" + id_selected + "')").eq(0).prop("selected", true);
            });
            $("#taxa").trigger("chosen:updated");
        }
    });
    return true;
}

function onRegionTypeChange(region_type){
    //var region_type = $("#select-region-type").val();
    if (region_type.value == "model"){
        $("#tr-model-select").show();
        // Finally, update gene models depending on v3 or v2
        //getGeneModels();  // TODO: move this line to onAssemblyChange();
        getRangeOfGeneModel($('#model-select').val(),$('#offset').val())
    }
    else{
        $("#tr-model-select").hide();
    }
}

function populateGeneModelTags(user_input){
    var matchingTags = [];
    if (user_input.length > 0){
        $.ajax({
            url: "get_gene_models.php",
            type: "POST",
            data: {
                action: "getGeneModels",
                assembly: $("#assembly").val(),
                input: user_input
            },
            success: function (data) {
                if (data == "null"){
                }
                else{
                    JSON.parse(data).forEach(function(obj) {
                        matchingTags.push(obj.model);
                    });
                    $( "#model-select" ).autocomplete({
                        source: matchingTags,
                        //minLength:2,
                        select: function(event, ui){
                            $( "#model-select" ).val(ui.item.value);
                        }
                    });

                    $( "#model-select" ).on("autocompleteselect", function(event, ui){
                        getRangeOfGeneModel(ui.item.value, $('#offset').val());
                    });
                }
            }
        });
    }
}

function getRangeOfGeneModel(gene_model, offset){
    if (gene_model != ""){
        $.ajax({
            type: "POST",
            data: {
                action: "getRange",
                assembly: $("#assembly").val(),
                model: gene_model
            },
            url: "get_gene_models.php",
            success: function (data) {
                var obj = JSON.parse(data);
                if (obj == null){
                    document.getElementsByClassName("verify-box")[0].innerHTML=
                        '<i class="fa fa-warning" aria-hidden="true" id="model-verify-box"></i>';
                }
                else{
                    document.getElementsByClassName("verify-box")[0].innerHTML=
                        '<i class="fa fa-check" aria-hidden="true" id="model-verify-box"></i>';
                    $("#chromosome").val(obj[0].chr);
                    $("#startPosition").val(Number(obj[0].min)-Number(offset));  // Start of gene model - offset
                    $("#endPosition").val(Number(obj[0].max)+Number(offset));  // End of gene model + offset
                    //$("#model-select").trigger("chosen:updated");
                }
                $(".verify-box").show();
            }
        });
    }
    else{
        $(".verify-box").hide();
    }
}

var config = {
    '.chosen-select': {},
    '.chosen-select-deselect': {allow_single_deselect: true},
    '.chosen-select-no-single': {disable_search_threshold: 10},
    '.chosen-select-no-results': {no_results_text: 'nothing found!'},
    '.chosen-select-width': {width: "95%"}
};
for (var selector in config) {
    $(selector).chosen(config[selector]);
}
