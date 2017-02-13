<html>
<body>
<?php 
$chromosome = $_POST["chromosome"];
$startPosition = $_POST["startPosition"]; 
$endPosition = $_POST["endPosition"]; 
$taxa = $_POST["taxa"];
$assembly = $_POST["dataSet"];
$version = $_POST["version"];

$Obj = array(
"chromosome" => $chromosome,
"startPosition" => $startPosition,
"endPosition" => $endPosition,
"taxa" => $taxa,
"dataSet" => $assembly,
"version" => $version,
);

$jsonObj=json_encode($Obj);
$$uId = "http://david1.usda.iastate.edu/Diversity/milestone2_service.php";
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
echo $response;
echo $status;

header('Location:'.$response);

echo "<a href=".$response.">Click here to view the output</a>";
?>
</body>
</html>
