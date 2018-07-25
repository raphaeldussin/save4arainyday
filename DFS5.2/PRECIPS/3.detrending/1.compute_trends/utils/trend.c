#include <stdio.h>
#include <string.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>
#include <netcdf.h>

/* Define the long names for netcdf variable */
#define trend_longname    "slope of the regression line"
#define correl_longname   "correlation time-variable"
#define autocor_longname  "autocorrelation (lag 1) of variable"
#define nfreedom_longname "effective degrees of freedom"
#define ttest_longname    "two-tailed t-test value"
#define fileout           "trend.nc"

#define NDIMSOUT 3

/* Handle errors by printing an error message and exiting with a
 *  * non-zero status. */
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

int
main(int argc, char *argv[])
{
   /* This will be the netCDF ID for the input file and data variable. */
   int ncid ; 
   int varid_data ;

   /* This will be the netCDF ID for the output file and data variable. */
   int ncid_out ;
   int x_dimid, y_dimid, timein_dimid ;
   int lonin_dimid, latin_dimid ;
   int timeout_dimid ;
   int dimids[NDIMSOUT]; 
   size_t start[NDIMSOUT], count[NDIMSOUT];
   int varid_trend, varid_correl, varid_autocor, varid_nfreedom, varid_ttest ;
   int varid_lonin, varid_latin, varid_lonout, varid_latout , varid_timeout ;
   /*char fileout[1024]; */

   /* Regression arrays and variables */
   double c0, c1, cov00, cov01, cov11, chisq;
   double r1, Neff, Nini, corr, test ;

   /* Loop indexes, and error handling. */
   int kx, ky, kt, retval;
   /*int dummy;*/
   int status=2;

   /* Strings from command line */
   int  iarg, ifile=0, nfiles=argc-13;
   char * variable = "", * time = "" ;
   char * dataset = "",  * fyear = ""; 
   char * lyear = "",    * diroutput = "" ; 
   char * listfile[nfiles] ;

   /* Arrays dimensions */
   size_t nframes, ntotal=0 ;
   size_t npx, npxsave=0, npy, npysave=0 ;
   char * clon="lon" , * clat="lat" ;
   int iframe ;

   /* time stuff */
   double zfyear, zlyear;
   double dummytime[1] ;

   dummytime[0] = 0. ;

   /* ------------------------------------------------------------- */
   /*                  Parse the command line                       */
   /* ------------------------------------------------------------- */

   for (iarg = 1; iarg < argc; iarg++){

       if (!strcmp(argv[iarg], "-var"))
          { variable = argv[iarg+1]; iarg = iarg +1; }
       else if (!strcmp(argv[iarg], "-time"))
          { time = argv[iarg+1]; iarg = iarg +1; }
       else if (!strcmp(argv[iarg], "-dataset"))
          { dataset = argv[iarg+1]; iarg = iarg +1; }
       else if (!strcmp(argv[iarg], "-fyear"))
          { fyear = argv[iarg+1]; iarg = iarg +1; }
       else if (!strcmp(argv[iarg], "-lyear"))
          { lyear = argv[iarg+1]; iarg = iarg +1; }
       else if (!strcmp(argv[iarg], "-diroutput"))
          { diroutput = argv[iarg+1]; iarg = iarg +1; }
       else
          { listfile[ifile] = argv[iarg]; ifile = ifile + 1 ; }
       }

   /* Verifying all fields are OK */
   if (strlen(variable) == 0)
      {printf ("Missing Variable... exiting \n")          ; exit(status) ; }
   if (strlen(time) == 0)
      {printf ("Missing time name... exiting \n")         ; exit(status) ; }
   if (strlen(dataset) == 0)
      {printf ("Missing Dataset... exiting \n")           ; exit(status) ; }
   if (strlen(fyear) == 0)
      {printf ("Missing First Year... exiting \n")        ; exit(status) ; }
   if (strlen(lyear) == 0)
      {printf ("Missing Last Year... exiting \n")         ; exit(status) ; }
   if (strlen(diroutput) == 0)
      {printf ("Missing Output Directory... exiting \n")  ; exit(status) ; }

   /* Control prints */
   printf(" Working on variable %s of dataset %s, from year %s to %s\n",variable, dataset, fyear, lyear);
   printf(" dimensions: time = %s \n",time);
   printf(" output directory = %s \n",diroutput);
   printf("\n");
   
   /*for (ifile=0; ifile < nfiles; ifile++){
       printf(" THE FILE %s\n", listfile[ifile]); } */

   /* ------------------------------------------------------------- */
   /*              Scan input files and allocate                    */
   /* ------------------------------------------------------------- */

   for (ifile=0; ifile < nfiles; ifile++){

      /* Open file */
      if ((retval = nc_open(listfile[ifile], NC_NOWRITE, &ncid)))
         ERR(retval);
      /* find time dimension id */
      if (( retval = nc_inq_dimid(ncid, time, &timein_dimid) ))
         ERR(retval);
      /* Find number of frames */
      if (( retval = nc_inq_dimlen(ncid, timein_dimid, &nframes) ))
         ERR(retval);

      /* Add to the total of frames */
      ntotal = ntotal + nframes ;

      /* find lon dimension id */
      if (( retval = nc_inq_dimid(ncid, clon, &lonin_dimid) ))
         ERR(retval);
      /* Find number of frames */
      if (( retval = nc_inq_dimlen(ncid, lonin_dimid, &npx) ))
         ERR(retval);
      /* Sanity check */
      if (( ifile == 0 ))
         npxsave = npx ;
      if (( npx != npxsave ))
         { printf(" lon dimension does not match between 2 files") ; exit(status) ; }

      /* find lat dimension id */
      if (( retval = nc_inq_dimid(ncid, clat, &latin_dimid) ))
         ERR(retval);
      /* Find number of frames */
      if (( retval = nc_inq_dimlen(ncid, latin_dimid, &npy) ))
         ERR(retval);
      /* Sanity check */
      if (( ifile == 0 ))
         npysave = npy ;
      if (( npy != npysave ))
         { printf(" lat dimension does not match between 2 files") ; exit(status) ; }

      /* Close the file */
      if ((retval = nc_close(ncid)))
         ERR(retval); }

   printf(" Nx = %ld , Ny = %ld , Nt = %ld",npxsave, npysave, ntotal);
   printf("\n");

   /* Define arrays size (waouh you can do that in C, surprising) */
   double data_in[ntotal][npysave][npxsave];

   double data_tmp[ntotal]; 
   double ztime[ntotal];

   double trend[1][npysave][npxsave], correl[1][npysave][npxsave];
   double autocor[1][npysave][npxsave], nfreedom[1][npysave][npxsave];
   double ttest[1][npysave][npxsave];
   double zlon[npxsave], zlat[npysave];

   /* ------------------------------------------------------------- */
   /*                   Read data from all files                    */
   /* ------------------------------------------------------------- */

   zfyear=strtod( fyear, NULL ); 
   zlyear=strtod( lyear, NULL ); 
   iframe=0;

   for (ifile=0; ifile < nfiles; ifile++){

      /* Open file */
      if ((retval = nc_open(listfile[ifile], NC_NOWRITE, &ncid)))
         ERR(retval);
      /* find time dimension id */
      if (( retval = nc_inq_dimid(ncid, time, &timein_dimid) ))
         ERR(retval);
      /* Find number of frames */
      if (( retval = nc_inq_dimlen(ncid, timein_dimid, &nframes) ))
         ERR(retval);

      double data_onefile[nframes][npysave][npxsave];

      if (( ifile == 0 ))
         { if ((retval = nc_inq_varid(ncid, clon, &varid_lonin)))
           ERR(retval); 
           if ((retval = nc_get_var_double(ncid, varid_lonin, &zlon[0])))
           ERR(retval);   
           if ((retval = nc_inq_varid(ncid, clat, &varid_latin)))
           ERR(retval);
           if ((retval = nc_get_var_double(ncid, varid_latin, &zlat[0])))
           ERR(retval);
           printf(" lon/lat read\n");
         }  
      /* Get the varid of the data variable, based on its name. */
      if ((retval = nc_inq_varid(ncid, variable, &varid_data)))
         ERR(retval);
      /* Read the data. */
      if ((retval = nc_get_var_double(ncid, varid_data, &data_onefile[0][0][0]))) 
         ERR(retval);
      /* Close the file, freeing all resources. */
      if ((retval = nc_close(ncid)))
         ERR(retval);

      /* Fill the big array */
      for (kt=0; kt<nframes; kt++){
          for (ky=0; ky<npysave; ky++){
              for (kx=0; kx<npxsave; kx++){
                  data_in[iframe+kt][ky][kx] = data_onefile[kt][ky][kx]; }}}

      iframe = iframe + nframes;

      }

      /* Fill time array */
      for (kt=0; kt<ntotal; kt++){
          ztime[kt] = zfyear + ( zlyear - zfyear + 1 ) * ( (double)kt / (double)ntotal ); }

      /*for (kt=0; kt<ntotal; kt++){
          printf("%f \n",ztime[kt]);} */

   /* ------------------------------------------------------------- */
   /*                            Main LOOP                          */
   /* ------------------------------------------------------------- */

   Nini = ntotal ;

   for (ky = 0; ky < npysave; ky++){
       for (kx = 0; kx < npxsave; kx++){
           for (kt = 0; kt < ntotal; kt++)
               { data_tmp[kt] = data_in[kt][ky][kx]; }

           /* Computing the linear trend for each point  */
           gsl_fit_linear (ztime, 1, data_tmp, 1, ntotal, 
                             &c0, &c1, &cov00, &cov01, &cov11, &chisq);

           /* Computing the correlation time/variable */
           corr = gsl_stats_correlation(ztime, 1, data_tmp, 1, ntotal) ;

           /* Computing the lag-1 autocorrelation */
           r1 = gsl_stats_lag1_autocorrelation(data_tmp, 1, ntotal) ;

           /* Original formula is : */
           /* r2 = 1.               */
           /* r1r2 = r1 * r2 ;      */
           /* Neff = Nini * (1-r1r2)/(1+r1r2) ; */
           /* But multiplication by one for each point is a waste of time */
           /* Effective degree of freedom */
           Neff = Nini * (1-r1)/(1+r1) ;

           /* compute the t-test */
           test = corr * sqrt((Neff-1)/(1-c1*c1)) ;

           /* Store the results into 2d-array for further writing */
           trend[0][ky][kx] = c1      ;
           correl[0][ky][kx] = corr   ;
           autocor[0][ky][kx] = r1    ;
           nfreedom[0][ky][kx] = Neff ;
           ttest[0][ky][kx] = test    ; }}

   /* ------------------------------------------------------------- */
   /*                    Create NetCDF output file                  */
   /* ------------------------------------------------------------- */

  /* dummy = sprintf(fileout, "%s/trend_%s_%s_%s-%s.nc",diroutput,variable,dataset,fyear,lyear); */
   /*printf("%s \n",fileout); */

   /* Create the file. The NC_CLOBBER parameter tells netCDF to
    * overwrite this file, if it already exists.*/
   if ((retval = nc_create(fileout, NC_CLOBBER, &ncid_out)))
      ERR(retval);

   /* Define the dimensions. NetCDF will hand back an ID for each. */
   if ((retval = nc_def_dim(ncid_out, "lon", npxsave, &x_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid_out, "lat", npysave, &y_dimid)))
      ERR(retval);
   if ((retval = nc_def_dim(ncid_out, "time", NC_UNLIMITED, &timeout_dimid)))
      ERR(retval);

   /* The dimids array is used to pass the IDs of the dimensions of
    * the variable. */
   dimids[0] = timeout_dimid; dimids[1] = y_dimid; dimids[2] = x_dimid;
   start[0] = 0; start[1] = 0; start[2] = 0;
   count[0] = 1; count[1] = npysave; count[2] = npxsave;

   /* Define the variables */ 
   if ((retval = nc_def_var(ncid_out, "lat", NC_DOUBLE, 1, &y_dimid, &varid_latout)))
      ERR(retval);

   if ((retval = nc_def_var(ncid_out, "lon", NC_DOUBLE, 1, &x_dimid, &varid_lonout)))
      ERR(retval);

   if ((retval = nc_def_var(ncid_out, "time", NC_DOUBLE, 1, &timeout_dimid, &varid_timeout)))
      ERR(retval);

   if ((retval = nc_def_var(ncid_out, "trend", NC_DOUBLE, NDIMSOUT, dimids, &varid_trend)))
      ERR(retval);

   if ((retval = nc_put_att_text(ncid_out, varid_trend, "long_name", strlen(trend_longname), trend_longname)))
      ERR(retval);

   if ((retval = nc_def_var(ncid_out, "correl", NC_DOUBLE, NDIMSOUT, dimids, &varid_correl)))
      ERR(retval);

   if ((retval = nc_put_att_text(ncid_out, varid_correl, "long_name", strlen(correl_longname), correl_longname)))
      ERR(retval);

   if ((retval = nc_def_var(ncid_out, "autocor", NC_DOUBLE, NDIMSOUT, dimids, &varid_autocor)))
      ERR(retval);

   if ((retval = nc_put_att_text(ncid_out, varid_autocor, "long_name", strlen(autocor_longname), autocor_longname)))
      ERR(retval);

   if ((retval = nc_def_var(ncid_out, "nfreedom", NC_DOUBLE, NDIMSOUT, dimids, &varid_nfreedom)))
      ERR(retval);

   if ((retval = nc_put_att_text(ncid_out, varid_nfreedom, "long_name", strlen(nfreedom_longname), nfreedom_longname)))
      ERR(retval);

   if ((retval = nc_def_var(ncid_out, "ttest", NC_DOUBLE, NDIMSOUT, dimids, &varid_ttest)))
      ERR(retval);

   if ((retval = nc_put_att_text(ncid_out, varid_ttest, "long_name", strlen(ttest_longname), ttest_longname)))
      ERR(retval);

   /* End define mode. */ 
   if ((retval = nc_enddef(ncid_out)))
      ERR(retval);

   /* Write coordinates */
   if (( retval = nc_put_var_double(ncid_out, varid_latout, &zlat[0]) ))
      ERR(retval);
   if (( retval = nc_put_var_double(ncid_out, varid_lonout, &zlon[0]) ))
      ERR(retval);
   if (( retval = nc_put_var_double(ncid_out, varid_timeout, &dummytime[0]) ))
      ERR(retval);

   /* Write the data to the file. */ 
   if ((retval = nc_put_vara_double(ncid_out, varid_trend, start, count,   &trend[0][0][0]))) 
      ERR(retval);
   if ((retval = nc_put_vara_double(ncid_out, varid_correl, start, count,  &correl[0][0][0]))) 
      ERR(retval);
   if ((retval = nc_put_vara_double(ncid_out, varid_autocor, start, count, &autocor[0][0][0]))) 
      ERR(retval);
   if ((retval = nc_put_vara_double(ncid_out, varid_nfreedom, start, count,&nfreedom[0][0][0]))) 
      ERR(retval);
   if ((retval = nc_put_vara_double(ncid_out, varid_ttest, start, count,   &ttest[0][0][0]))) 
      ERR(retval);

   /* Close the file. */
   if ((retval = nc_close(ncid_out)))
      ERR(retval);

   return 0;
}
