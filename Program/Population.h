/*  ---------------------------------------------------------------------- //
    Hybrid Genetic Search for Arc Routing Problems -- HGS-CARP
    Copyright (C) 2016 Thibaut VIDAL

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//  ---------------------------------------------------------------------- */

#ifndef POPULATION_H
#define POPULATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include "Noeud.h"
#include "Individual.h"

using namespace std ;

// Structure to manage a sub-population (feasible or infeasible solutions)
struct SubPop {
    // Individuals
    vector<Individual *> individus;

    // Number of Individuals
    int nbIndiv;
};

class Population
{
   private:

   // List to remember which of the 50 previous individuals were feasible in terms of load capacity
   list <bool> listeValiditeCharge ;

   // List to remember which of the 50 previous individuals were feasible in terms of max travel time
   list <bool> listeValiditeTemps ;

   // Education procedure (LS)
   void education(Individual * indiv);

    // Place an individual in the population
    // Returns its position
    int placeIndividu(SubPop *pop, Individual *indiv);

   public:

   // Access to the parameters of the problem
   Params * params ;

   // clock time when the best individual was found
   clock_t timeBest ;

   // Auxiliary data structure (Individual) with all local search data structures
   // To do the LS on a given individual, we simply copy in this individual and run the LS there.
   Individual * trainer;

    // check if there is already a solution with the same fitness in the population
    bool fitExist(SubPop *pop, Individual *indiv);

    // compute the biased fitness of the individuals in the population
    void evalExtFit(SubPop *pop);

   // add an individual in the population
   int addIndividu (Individual * indiv);

   // add all individuals from another population
   int addAllIndividus (Population * pop);

    // remove an individual in the population (chosen accoding to the biased fitness)
    void removeIndividu(SubPop *pop, int p);

    // subprocedure that chooses an individual to be removed
    int selectCompromis(SubPop *souspop);

    // update the table of distance (Hamming distance) between individuals to know their proximity
    void updateProximity(SubPop *pop, Individual *indiv);

    // Diversification procedure (replace a large part of the population by new random solutions)
    void diversify();

    // Clear (no more individual in the population), used only in the ILS version of the code
    void clear();

    // Feasible and Infeasible subpopulations
    SubPop *valides;
    SubPop *invalides;

    // Get one individual per binary tournament
    Individual *getIndividuBinT();

    // Get one individual with uniform probability in a percentage of the best
    Individual *getIndividuPourc(int pourcentage);

    // Get best feasible individual
    Individual *getIndividuBestValide();

   // Get best infeasible individual
   Individual * getIndividuBestInvalide ();

    // when the penalty coefficient change, need to recompute properly the fitness of the individuals in the population
    void validatePen(SubPop *subPop);

   //////////////////////////////////////////////////////////

   // Print the best solution in a file
   void ExportBest (string filename) ;

   // Solution check
   // Verifies the cost and feasibility of the solution
   // Using only the instance data and the shortest path data (not relying on the auxiliary data structures)
   bool solutionChecker(vector < vector < vector < int > > > & allRoutes, vector < vector < vector < pair <int,int > > > > & allRoutesArcs, double expectedCost, double expectedMaxRoute);
   
   // Print the best solution in the BKS file, only if its better than the previous BKS
   void ExportBKS (string nomFichier) ;

   /* FUNCTIONS FOR REGULAR PRINTOUTS OF THE STATUS OF THE POPULATION */

   // get the fraction of valid individuals with respect to the load constraint
   double fractionValidesCharge () ;

   // get the fraction of valid individuals with respect to the time constraint
   double fractionValidesTemps () ;

    // get the diversity of the population
    double getDiversity(SubPop *pop);

   // get the average cost of a feasible solution in the population
   double getMoyenneValides ();

   // get the average cost of an infeasible solution in the population
   double getMoyenneInvalides ();

   // get the average age of feasible solutions
   double getAgeValides ();

   // print a small report of the status of the population
   void afficheEtat(int NbIter);

   // update the count of valid individuals
   void updateNbValides (Individual * indiv);

   // update the age of the individuals
   void updateAge ();

   // Constructor
   Population(Params * params);

   // Destructor
   ~Population();
};

#endif
