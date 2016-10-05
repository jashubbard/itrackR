#include <stdio.h>
#include <sqlite3.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// callback function when running SELECT query. Prints to screen
static int callback(void *NotUsed, int argc, char **argv, char **azColName){
  int i;
  for(i=0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  printf("\n");
  return 0;
}


// [[Rcpp::export]]
int fast_sql(std::string dbname, std::string sql)
{
  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;

  // convert from std::string to const char * for sqlite3_open
  const char * dbn = dbname.c_str();

  rc = sqlite3_open(dbn, &db);

  if( rc ){
    printf("Can't open database: %s\n", sqlite3_errmsg(db));
    return(0);
  }else{
    fprintf(stderr, "Opened database successfully\n");
  }

  // convert to const char for sqlite3_exec
  const char *query = sql.c_str();

  /* Execute SQL statement */
  rc = sqlite3_exec(db, query, callback, 0, &zErrMsg);

  // give error/success message
  if( rc != SQLITE_OK ){
    Rprintf("SQL error: %s\n", zErrMsg);
    sqlite3_free(zErrMsg);
  }else{
    // fprintf(stdout, "Records created successfully\n");
  }

  // close the db
  sqlite3_close(db);

  return(rc);
}

// [[Rcpp::export]]
int update_rows(std::string dbname, NumericVector tstart, NumericVector tend){

  const char *fmt;
  char *sql = (char*)malloc(13 * sizeof(char));
  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;
  const char *query;

  // convert from std::string to const char * for sqlite3_open
  const char * dbn = dbname.c_str();

  rc = sqlite3_open(dbn, &db);

  if( rc ){
    printf("Can't open database: %s\n", sqlite3_errmsg(db));
    return(0);
  }else{
    fprintf(stderr, "Opened database successfully\n");
  }

  rc = sqlite3_exec(db, "PRAGMA synchronous=OFF; PRAGMA journal_mode=MEMORY", callback, 0, &zErrMsg);


  for(int i=0; i<tstart.size(); i++){

  /* Create merged SQL statement */
  // fmt = "UPDATE SAMPLES set TEST = %d where time>=%d " \
  // "AND time <= %d; ";

  fmt = "SELECT * FROM SAMPLES where time>=%d " \
    "AND time <= %d; ";

  sprintf(sql,fmt,i,tstart[i],tend[i]);


  /* Execute SQL statement */
  rc = sqlite3_exec(db, sql, callback, 0, &zErrMsg);

  // give error/success message
  if( rc != SQLITE_OK ){
    Rprintf("SQL error: %s\n", zErrMsg);
    sqlite3_free(zErrMsg);
  }

  }

  // close the db
  sqlite3_close(db);

  return(0);

}

//
// arma::mat update_df(NumericMatrix samps, NumericVector starttimes, NumericVector endtimes,NumericVector vals){
//
//   arma::vec st=arma::vec(starttimes);
//   arma::vec en=arma::vec(endtimes);
//   arma::vec val=arma::vec(vals);
//
//
//
//
//
// }


