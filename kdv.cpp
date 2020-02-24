#include <iostream>
#include <iomanip>
#include <math.h>
#include <getopt.h>
#include "kdvintegrator.h"

// Default parameter values
const double      C_DEFAULT_XMAX = 50.0;
const double      C_DEFAULT_TMAX = 10.0;
const double      C_DEFAULT_DX   = 0.1;
const double      C_DEFAULT_DT   = 0.00001;
const std::string C_DEFAULT_FILE = "kdv_out";


static void usage(int argc, char *argv[])
{
   if (argc>0)
   {
      fprintf(stderr, "Usage: %s\n", argv[0]);
   }
   fprintf(stderr, "\t--xmax <xmax>\t(default: %f)\n", C_DEFAULT_XMAX);
   fprintf(stderr, "\t--tmax <tmax>\t(default: %f)\n", C_DEFAULT_TMAX);
   fprintf(stderr, "\t--dx   <dx>  \t(default: %f)\n", C_DEFAULT_DX);
   fprintf(stderr, "\t--dt   <dt>  \t(default: %f)\n", C_DEFAULT_DT);
   fprintf(stderr, "\t--file <file>\t(default: %s)\n", C_DEFAULT_FILE.c_str());
}


int main(int argc, char *argv[])
{
   double xmax      = C_DEFAULT_XMAX;
   double tmax      = C_DEFAULT_TMAX;
   double dx        = C_DEFAULT_DX;
   double dt        = C_DEFAULT_DT;
   std::string file = C_DEFAULT_FILE;

   int c;
   while (1)
   {
      int option_index = 0;
      static struct option long_options[] =
      {
         {"xmax", required_argument, 0, 0},
         {"tmax", required_argument, 0, 0},
         {"dx",   required_argument, 0, 0},
         {"dt",   required_argument, 0, 0},
         {"file", required_argument, 0, 0},
         {0,      0,                 0, 0}
      };

      c = getopt_long(argc, argv, "", long_options, &option_index);
      if (c == -1)
         break;

      if (c == 0)
      {
         switch (option_index)
         {
            case 0 : xmax = atof(optarg); break;
            case 1 : tmax = atof(optarg); break;
            case 2 : dx   = atof(optarg); break;
            case 3 : dt   = atof(optarg); break;
            case 4 : file = optarg; break;
            default: usage(argc, argv);
                     exit(EXIT_FAILURE);
         }
      } else {
         usage(argc, argv);
         exit(EXIT_FAILURE);
      }
   } // end of while (1)

   std::cout << "Running KdV simulation with the following parameters:" << std::endl;
   std::cout << "xmax : " << std::setw(19) << std::setprecision(4) << xmax << std::endl;
   std::cout << "tmax : " << std::setw(19) << std::setprecision(4) << tmax << std::endl;
   std::cout << "dx   : " << std::setw(19) << std::setprecision(4) << dx << std::endl;
   std::cout << "dt   : " << std::setw(19) << std::setprecision(4) << dt << std::endl;
   std::cout << "file : " << file << std::endl;


   const unsigned int nmax = xmax/dx + 0.5;  // Add 0.5 for rounding.
   const unsigned int mmax = tmax/dt + 0.5;  // Add 0.5 for rounding.

   dv u(nmax+1);

   KdvIntegrator kdvintegrator(mmax, nmax, dt, dx, file);

   kdvintegrator.initialize(u);
   for (unsigned int m=0; m<=mmax; ++m)
   {
      const double t = m*dt;
      kdvintegrator.single_step(t, u);
   }

   return 0;
} // main

