<!-- Web client for second use case, with smaller list of lines/taxa. 
Sends input via POST-method to milestone2_send.php
Use-case: "Given a range of genomic positions and a specific line, what sequence of allele is present?
          List all the lines which have the same sequence in the same range of genomic positions."-->
<html>
	<head>
		<title> Query Genotype Data</title>
		<script type="text/javascript" src="js/jquery.js"></script>
		<script type="text/javascript" src="js/home.js"></script>
	</head>
	<body>
	 	<form action="milestone2_send.php" method="post">
		<table border="2">
			<tr>
				<td> Chromosome : 
				</td>
				<td>
					<select id="chromosome" name="chromosome">
					  <option value="0">0</option>
					  <option value="1">1</option>
					  <option value="2">2</option>
					  <option value="3">3</option>
					  <option value="4">4</option>
					  <option value="5">5</option>
					  <option value="6">6</option>
					  <option value="7">7</option>
					  <option value="8">8</option>
					  <option value="9">9</option>
					  <option value="10">10</option>
					</select>
				</td>
			</tr>
			<tr>
				<td> Positions :
				</td>
				<td> Start : <input id="startPosition" name="startPosition" type="number"> 
					End : <input id="endPosition" name="endPosition" type="number">
				</td>
			</tr>
			<tr>
				<td>
				Select Taxa :
				</td>
				<td>
					 <select id="taxa" name="taxa">
					 	<option value='Ames10244:250032849'>Ames10244:250032849</option>
						<option value='Ames10246:250041697'>Ames10246:250041697</option>
						<option value='Ames10247:250033161'>Ames10247:250033161</option>
						<option value='Ames10248:250031797'>Ames10248:250031797</option>
						<option value='Ames10249:250033148'>Ames10249:250033148</option>
						<option value='Ames10250:250033048'>Ames10250:250033048</option>
						<option value='Ames10251:250033046'>Ames10251:250033046</option>
						<option value='Ames10252:250033235'>Ames10252:250033235</option>
						<option value='Ames10253:250031774'>Ames10253:250031774</option>
						<option value='Ames10254:250033198'>Ames10254:250033198</option>
						<option value='Ames10255:250033271'>Ames10255:250033271</option>
						<option value='Ames10256:250033096'>Ames10256:250033096</option>
						<option value='Ames10256:250042417'>Ames10256:250042417</option>
						<option value='Ames10257:250033171'>Ames10257:250033171</option>
						<option value='Ames10258:250033159'>Ames10258:250033159</option>
						<option value='Ames10259:250033213'>Ames10259:250033213</option>
						<option value='Ames10260:250031197'>Ames10260:250031197</option>
						<option value='Ames10261:250033107'>Ames10261:250033107</option>
						<option value='Ames10262:250033038'>Ames10262:250033038</option>
						<option value='Ames10263:250031024'>Ames10263:250031024</option>
						<option value='Ames10264:250031258'>Ames10264:250031258</option>
						<option value='Ames10265:250035124'>Ames10265:250035124</option>
						<option value='Ames10266:250042451'>Ames10266:250042451</option>
						<option value='Ames10267:250031659'>Ames10267:250031659</option>
						<option value='Ames10268:250042490'>Ames10268:250042490</option>
						<option value='Ames10269:250033220'>Ames10269:250033220</option>
						<option value='Ames10271:250031710'>Ames10271:250031710</option>
						<option value='Ames10272:250033211'>Ames10272:250033211</option>
						<option value='Ames10273:250031771'>Ames10273:250031771</option>
						<option value='Ames10274:250033192'>Ames10274:250033192</option>
						<option value='Ames10274:250042502'>Ames10274:250042502</option>
						<option value='Ames10276:250033172'>Ames10276:250033172</option>
						<option value='Ames10277:250033095'>Ames10277:250033095</option>
						<option value='Ames10287:250032568'>Ames10287:250032568</option>
						<option value='Ames12725:250034925'>Ames12725:250034925</option>
						<option value='Ames12726:250032944'>Ames12726:250032944</option>
						<option value='Ames12727:250035100'>Ames12727:250035100</option>
						<option value='Ames12728:250032927'>Ames12728:250032927</option>
						<option value='Ames12729:250033110'>Ames12729:250033110</option>
						<option value='Ames12729:250033153'>Ames12729:250033153</option>
						<option value='Ames12730:250033034'>Ames12730:250033034</option>
						<option value='Ames12731:250032800'>Ames12731:250032800</option>
						<option value='Ames12732:250032896'>Ames12732:250032896</option>
						<option value='Ames12732:250032907'>Ames12732:250032907</option>
						<option value='Ames12733:250032773'>Ames12733:250032773</option>
						<option value='Ames12734:250032780'>Ames12734:250032780</option>
						<option value='Ames12734:250032850'>Ames12734:250032850</option>
						<option value='Ames12735:250033219'>Ames12735:250033219</option>
						<option value='Ames12736:250031951'>Ames12736:250031951</option>
						<option value='Ames12737:250032770'>Ames12737:250032770</option>
						<option value='Ames12738:250032926'>Ames12738:250032926</option>
						<option value='Ames12739:250033075'>Ames12739:250033075</option>
						<option value='Ames12740:250033043'>Ames12740:250033043</option>
						<option value='Ames12814:250031891'>Ames12814:250031891</option>
						<option value='Ames12815:250031896'>Ames12815:250031896</option>
						<option value='Ames12816:250041743'>Ames12816:250041743</option>
						<option value='Ames12817:250031903'>Ames12817:250031903</option>
						<option value='Ames12817:250041694'>Ames12817:250041694</option>
						<option value='Ames12818:250033615'>Ames12818:250033615</option>
						<option value='Ames12818:250041678'>Ames12818:250041678</option>
						<option value='Ames12819:250031873'>Ames12819:250031873</option>
						<option value='Ames12820:250033583'>Ames12820:250033583</option>
						<option value='Ames12821:250031934'>Ames12821:250031934</option>
						<option value='Ames12822:250033599'>Ames12822:250033599</option>
						<option value='Ames12823:250033553'>Ames12823:250033553</option>
						<option value='Ames12824:250031960'>Ames12824:250031960</option>
						<option value='Ames14111:250033032'>Ames14111:250033032</option>
						<option value='Ames14111:250041653'>Ames14111:250041653</option>
						<option value='Ames14112:250031426'>Ames14112:250031426</option>
						<option value='Ames14113:250032833'>Ames14113:250032833</option>
						<option value='Ames14113:250042374'>Ames14113:250042374</option>
						<option value='Ames14114:250032791'>Ames14114:250032791</option>
						<option value='Ames14115:250034754'>Ames14115:250034754</option>
						<option value='Ames14116:250033794'>Ames14116:250033794</option>
						<option value='Ames14211:250031517'>Ames14211:250031517</option>
						<option value='Ames14408:250031965'>Ames14408:250031965</option>
						<option value='Ames15929:250033855'>Ames15929:250033855</option>
						<option value='Ames18999:250032813'>Ames18999:250032813</option>
						<option value='Ames18999:250032902'>Ames18999:250032902</option>
						<option value='Ames19000:250031259'>Ames19000:250031259</option>
						<option value='Ames19000:250031373'>Ames19000:250031373</option>
						<option value='Ames19000:250041598'>Ames19000:250041598</option>
						<option value='Ames19000:250042561'>Ames19000:250042561</option>
						<option value='Ames19001:250034888'>Ames19001:250034888</option>
						<option value='Ames19002:250033017'>Ames19002:250033017</option>
						<option value='Ames19002:250033342'>Ames19002:250033342</option>
						<option value='Ames19003:250033145'>Ames19003:250033145</option>
						<option value='Ames19004:250033179'>Ames19004:250033179</option>
						<option value='Ames19005:250032014'>Ames19005:250032014</option>
						<option value='Ames19006:250031408'>Ames19006:250031408</option>
						<option value='Ames19007:250031290'>Ames19007:250031290</option>
						<option value='Ames19008:250031330'>Ames19008:250031330</option>
						<option value='Ames19008:250041562'>Ames19008:250041562</option>
						<option value='Ames19009:250033582'>Ames19009:250033582</option>
						<option value='Ames19010:250034541'>Ames19010:250034541</option>
						<option value='Ames19011:250031433'>Ames19011:250031433</option>
						<option value='Ames19012:250033821'>Ames19012:250033821</option>
						<option value='Ames19013:250033793'>Ames19013:250033793</option>
						<option value='Ames19014:250033702'>Ames19014:250033702</option>
						<option value='Ames19015:250031558'>Ames19015:250031558</option>
					 </select>
				</td>
			<tr>
				<td> Data Set:
				</td>
				<td>
					<select id="dataSet" name="dataSet">
					  <option value="amesTest">Ames Test</option>
					  <option value="AllZeaGBSv27public20140528">AllZeaGBSv27public20140528</option>
					  <option value="ZeaGBSv27publicImputed20150114">ZeaGBSv27publicImputed20150114</option>
					</select>
				</td>
			</tr>
				<td> Version:
				</td>
				<td>
					<select id="version" name="version">
					  <option value="stable">stable</option>
					  <option value="alpha">alpha</option>
					</select>
				</td>
			</tr>
			<tr>
				<td colspan="2">
					<input type="submit">
				</td>
			</tr>
			
		</table>
		</form>
	<img id="loading" src="img/ajax-loader.gif" style="display: none">
	</body>
</html>
