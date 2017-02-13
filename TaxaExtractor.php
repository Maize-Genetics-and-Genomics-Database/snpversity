<?php

require_once 'db_controller.php';
/* Used for extracting and preprocessing all taxa to ensure appropriate input for tassel.
 * data is retrieved from send.php, processed, and array of taxa is send back.
 */

class TaxaExtractor {

    private $taxaArray;
    private $db_handle;
    private $projsAsTaxon;
    private $taxaFileName;
    private $taxaFileTemp;
    private $projects;
    private $assembly;

    function __construct($assembly) {
        $this->taxaArray = array();
        $this->projsAsTaxon = array();
        $this->projects = array();
        $this->db_handle = new DBController();
        $this->assembly = $assembly;
    }

    public function setTaxaArray($taxa) {
        $this->taxaArray = $taxa;
    }

    public function setTaxaDataWithFiles($taxa, $taxaFileName, $taxaFile) {
        $this->taxaArray = $taxa;
        $this->taxaFileName = $taxaFileName;
        $this->taxaFileTemp = $taxaFile;
    }

    public function extract() {
        //global $db_handle, $projsAsTaxon,$taxaFileName, $taxaFileTemp, $taxaArray, $projects;
        $query = "SELECT DISTINCT project FROM allzeagbsv27;";
        $results = $this->db_handle->runQuery($query);

        //Load project names into $projects
        foreach ($results as $row) {
            array_push($this->projects, $row["project"]);
        }
        //could have submitted empty taxaArray
        if (!isset($this->taxaArray)) {
            $this->taxaArray = array();
        }
        //Find out if user selected 'all taxa in project' (taxon-name = project)
        $this->projsAsTaxon = array_intersect($this->projects, $this->taxaArray);

        //projects could be in file
        if (!isset($this->projsAsTaxon)) {
            $this->projsAsTaxon = array();
        }
        $extension = $this->getFileExtension($this->taxaFileName);
        if ($extension == "stockinfo" ||$extension == "txt" || $extension == "taxainfo") {
            $this->getTaxaFromTaxaInfoFile($this->taxaFileTemp);
        } elseif ($extension == "csv") {
            $this->addTaxaFromCsvFile($this->taxaFileTemp);
        }
        $this->addAllTaxaFromProj($this->projsAsTaxon);
        return array_filter($this->taxaArray);
    }

    private function addAllTaxaFromProj($selectedProjs) {
        foreach ($selectedProjs as $proj) {
            if ($proj == "RandD") {
                $proj = "R&D";
            }
            if ($proj == "All") {
                $query = "SELECT DISTINCT ON (dna_sample,project) dna_sample, lib_prep_id "
                        . "FROM allzeagbsv27";
            } else {
                $query = "SELECT DISTINCT ON (dna_sample) dna_sample, lib_prep_id "
                        . "FROM allzeagbsv27 "
                        . "WHERE project = '" . $proj . "' "
                        . "ORDER BY dna_sample, lib_prep_id";
            }
            $results = $this->db_handle->runQuery($query);
            foreach ($results as $row) {
                $taxon = $row["dna_sample"] . ":" . $row["lib_prep_id"];
                if (!in_array($taxon, $this->taxaArray) && (strlen($taxon) >0)) { // handle duplicates
                    array_push($this->taxaArray, $taxon);
                }
            }
            //Finally, remove project from taxa-selections
            if (($key = array_search($proj, $this->taxaArray)) !== false) {
                unset($this->taxaArray[$key]);
            }
        }
    }

    public function getTaxaFromTaxaInfoFile($fileName) {
        // Add taxa from file, if submitted
        if (!empty($fileName)) {
            $taxaArrayFromFile = file($fileName, FILE_IGNORE_NEW_LINES);
            foreach ($taxaArrayFromFile as $taxonToAdd) {
                $this->extractTaxaFromLine($taxonToAdd);
            }
            if (end($this->taxaArray) == "") { //dont want newline at the end
                array_pop($this->taxaArray);
            }
        }
        return array_filter($this->taxaArray);
    }

    public static function getFileExtension($filename) {
        if (!empty($filename)) {
            $parts = explode('.', $filename);
            return $parts[count($parts) - 1];
        } else {
            return "";
        }
    }

//extracts any taxa from CSV file and pushes them to global $taxaArray
    private function addTaxaFromCsvFile($fileName) {
        // Add taxa from file, if submitted
        if (!empty($fileName)) {
            $fp = fopen($fileName, "r");
            if ($fp !== false) {
                while (!feof($fp)) {
                    $row = fgetcsv($fp, 1024);
                    $cols_len = count($row);
                    if ($cols_len == 1) {
                        $this->extractTaxaFromLine($row[0]);  //treat like .taxainfo file because only one col
                    } elseif ($cols_len == 2) {  //have project, taxon
                        $proj = $row[0];
                        $taxon_commented = preg_replace('/\s+/', '', $row[1]);
                        //drop part after #, because comment
                        $taxon_split = explode('#', $taxon_commented);
                        $taxon = $taxon_split[0];
                        $taxon = $this->checkIfInbred($taxon,true); // check if inbred
                        if ($proj == "All" || $proj == "all") { //look through all projects to find match
                            $query = "SELECT DISTINCT ON (project) dna_sample, lib_prep_id"
                                    . " FROM allzeagbsv27 "
                                    . "WHERE dna_sample = '" . $taxon . "' ORDER BY project, lib_prep_id;";
                            $results = $this->db_handle->runQuery($query);
                            $this->addToTaxaArray($results, false);
                        } elseif (($taxon == "All" || $taxon == "all") && in_array($proj, $this->projects)) { //Case Project, All <-> select all taxa from a project
                            if (!in_array($proj, $this->projsAsTaxon)) {
                                array_push($this->projsAsTaxon, $proj);
                            }
                        } else {  //case project,taxon
                            $query = "SELECT DISTINCT ON (dna_sample) dna_sample, lib_prep_id"
                                    . " FROM allzeagbsv27 "
                                    . "WHERE dna_sample = '" . $taxon . "' AND project='" . $proj . "' ORDER BY dna_sample, lib_prep_id;";
                            $results = $this->db_handle->runQuery($query);
                            $this->addToTaxaArray($results, false);
                        }
                    }
                }
            }
        }
    }

    private function addToTaxaArray($results, $taxonIsHapMap) {
        if ($this->assembly == "v2") { // !HapMap
            foreach ($results as $row) {
                $taxonToAddFiltered = $row["dna_sample"] . ':' . $row["lib_prep_id"];
                if (!in_array($taxonToAddFiltered, $this->taxaArray) && !(empty($taxonToAddFiltered))) {
                    array_push($this->taxaArray, $taxonToAddFiltered);
                }
            }
        } else {
            foreach ($results as $row) {
                $taxon = $row["taxon"];
                if (!in_array($taxon, $this->taxaArray) && !empty($taxon)) {
                    array_push($this->taxaArray, $taxon);
                }
            }
        }
    }

    //extract all taxa from one line in .stockinfo file
    private function extractTaxaFromLine($taxonToAdd) {
        $taxonToAddFiltered = str_replace(array("\r", "\n"), "", $taxonToAdd);  //remove newlines and empty spaces
        $taxonToAddFiltered = $this->checkIfInbred($taxonToAddFiltered);
        //$taxonNotAllZeaGBS = (substr_count($taxonToAddFiltered, ':') == 0); // does not have a ':' and is not a project (only HapMap dataset)

        //handle common name
        if ($this->assembly == "v3") {//!= allzeagbsv27 dataset
            //get max 1 result per project
            $results = $this->db_handle->runQuery("SELECT * "
                    . "FROM hapmapv3 "
                    . "WHERE taxon = '" . $taxonToAddFiltered . "' "
                    . "ORDER BY taxon;");
            $taxonToAddFiltered = "";
            $this->addToTaxaArray($results, true);
        }

        // treat as dna_sample:lib_prep_id (same as generated file)
        elseif (!in_array($taxonToAddFiltered, $this->taxaArray) && !empty($taxonToAddFiltered)) {
            array_push($this->taxaArray, $taxonToAddFiltered);
        }
    }

//extract all taxa from one line in .taxainfo file [DEPRECATED]
   /* private function extractTaxaFromLineOld($taxonToAdd, $taxaArray) {
        $taxonToAddFiltered = str_replace(array("\r", "\n"), "", $taxonToAdd);  //remove newlines and empty spaces
        $taxonIsProj = in_array($taxonToAddFiltered, $this->projects);
        $taxonFromHapMap = (substr_count($taxonToAddFiltered, ':') == 0) && !$taxonIsProj; // does not have a ':' and is not a project (only HapMap dataset)
        /* handle full names [deprecated]
          if (substr_count($taxonToAddFiltered, ':') == 3) {  //full name
          $taxonToAddFiltered = (preg_replace("/:(.*?):\d/", "", $taxonToAddFiltered));
          } */
        /*handle common name
        if ($taxonFromHapMap) {  //get max 1 result per project
            $results = $this->db_handle->runQuery("SELECT * "
                    . "FROM hapmapv3 "
                    . "WHERE taxon = '" . $taxonToAddFiltered . "' "
                    . "ORDER BY taxon;");
            $taxonToAddFiltered = "";
            $this->addToTaxaArray($results, $taxonFromHapMap);
        }
        // treat as dna_sample:lib_prep_id (same as generated file)
        elseif (!in_array($taxonToAddFiltered, $taxaArray) && !$taxonIsProj) {
            array_push($taxaArray, $taxonToAddFiltered);
        } elseif ($taxonIsProj) {
            if (!in_array($taxonToAddFiltered, $this->projsAsTaxon)) {
                array_push($this->projsAsTaxon, $taxonToAddFiltered);
            }
        }
    }*/

    public function checkIfInbred($inbred, $shortName=false){
        if ($this->assembly == "v2"){  // Currently only have inbreds for allzeagbs
            $query_get_accs = "SELECT DISTINCT ON (dna_sample,project) * FROM allzeagbsv27 WHERE inbred='".$inbred."';";
            $accessions = $this->db_handle->runQuery($query_get_accs);
            if ($accessions[0]["dna_sample"] != ''){
                if ($shortName){
                    return $accessions[0]["dna_sample"];
                }
                return $accessions[0]["dna_sample"] . ":" .$accessions[0]["lib_prep_id"];
            }
        }
        return $inbred;  // nothing found
    }
}
