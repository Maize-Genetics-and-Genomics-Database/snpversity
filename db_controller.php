<?php
class DBController {
        private $host = "host=127.0.0.1";
        private $port = "port=5432";
        private $dbname = "dbname=postgres";
        private $credentials = "user=postgres password=david101";  //TODO: Encrypt password and user?
	function __construct() {
		$conn = $this->connectDB();
	}
	
	function connectDB() {
		$conn = pg_connect($this->host . " ".$this->port. " " . $this->dbname . " " . $this->credentials);
		return $conn;
	}
	
	function runQuery($query) {
		$result = pg_query($query);
		while($row=pg_fetch_assoc($result)) {
			$resultset[] = $row;
		}		
		if(!empty($resultset))
			return $resultset;
	}
	
	function numRows($query) {
		$result  = pg_query($query);
		$rowcount = pg_num_rows($result);
		return $rowcount;	
	}
}
?>