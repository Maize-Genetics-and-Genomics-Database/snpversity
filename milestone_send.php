<html>
<body>
<?php 
$chromosome = $_POST["chromosome"];
$position = $_POST["position"]; 
$taxa = $_POST["taxa"];
$assembly = $_POST["dataSet"];
$version = $_POST["version"];

$Obj = array(
"chromosome" => $chromosome,
"position" => $position,
"taxa" => $taxa,
"dataSet" => $assembly,
"version" => $version,
);

$jsonObj=json_encode($Obj);
$$uId = "http://david1.usda.iastate.edu/Diversity/milestone_service.php";
$curl=curl_init();

curl_setopt_array($curl, array(
    CURLOPT_RETURNTRANSFER => 1,
    CURLOPT_URL => $$uId,
    CURLOPT_USERAGENT => 'Codular Sample cURL Request',
    CURLOPT_POST => 1,
    CURLOPT_POSTFIELDS => array(
        data => $jsonObj
    )
));

$response = curl_exec($curl);
$status = curl_getinfo($curl, CURLINFO_HTTP_CODE);
curl_close($curl);

header('Location:'.$response);

echo "<a href=".$response.">Click here to view the output</a>";
?>
</body>
</html>
