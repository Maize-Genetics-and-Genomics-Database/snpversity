$(function () {
    $("#fileModal").load("html/file_modal.html");
    $("#helpModal").load("html/help_modal.html");
    $("#hapmapV3Modal").load("html/hapmapV3_modal.html");
    $("#taxaSelectModal").load("html/taxa_select_modal.html");
    $("#allzeagbsModal").load("html/allzeagbs_modal.html");
    $("#geneModelsModal").load("html/gene_models_modal.html");
    $("#datasetsModal").load("html/datasets_modal.html");
});

//deprecated
function addHelp(id) {
    var html = "<div id = " + id + "></div>";
    var data = $(id).load("html/" + id + ".html");
    html.innerHTML = data;
    return html;
}