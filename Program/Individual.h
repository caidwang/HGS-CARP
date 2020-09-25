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

#ifndef INDIVIDU_H
#define INDIVIDU_H

#include <vector>
#include <list>
#include <iostream> 
#include "Noeud.h"
#include "Params.h"
#include "LocalSearch.h"
using namespace std ;
class LocalSearch ;

// Structure used to keep the cost of a solution
struct CostSol {

	// value of the fitness, including the penalties for possible violations of the constraints
	double evaluation ;

	// value of the distance (not including the penalty values)
	double distance ;

	// amount of capacity violation
	double capacityViol ;

	// amout of distance violation
	double lengthViol ;

	// number of routes
	int routes ;


	// Constructor
	CostSol(){
		evaluation = 0 ;
		distance = 0 ;
		capacityViol = 0 ;
		lengthViol = 0 ;
		routes = 0 ;
	}

	// says if "this" is better in terms of feasibility than "c2"
	bool isBetterFeas(CostSol c2){
		bool isBetter = false ;
		bool isDominated = false ;
		if (capacityViol > c2.capacityViol + 0.0001) isDominated = true ;
		if (lengthViol > c2.lengthViol + 0.0001) isDominated = true ;
		if (capacityViol < c2.capacityViol - 0.0001) isBetter = true ;
		if (lengthViol < c2.lengthViol - 0.0001) isBetter = true ;
		return (isBetter && !isDominated);
	}
};

// preliminary declaration, because proxData depends on Individual
class Individual ;

// data about the proximity of an Individual with regards to the others in the population.
struct proxData {
	// Individual
	Individual * indiv ;

	// Its distance to the others
	double dist ;
};

class Individual
{

public:

	// Access to the parameters of the problem
	Params * params ;

	// How old is this individual
	int age ;

	// fitness of the individual
	double fitnessExtended ;

	// rank in terms of diversity
	float divRank ;

	// rank in terms of cost
	float fitRank ;

	// evaluation of the solution
	CostSol costSol ;

	// number of routes (if computed)
	int nbRoutes ;

	// length of the maximum route (if computed, the computation is done during split)
	double maxRoute ;

	// The giant tour of the individual
	// chromT[i][j] -> day i, j^th customer in the giant tour
	// 保存giant tour的内容
	vector < vector<int> > chromT ;

	// Pattern chromosome of an individual
	// 保存pattern的内容和分配的depot
	vector < pattern > chromP ;

	// Indices of the beginning of the routes (if already computed via Split)
	// chromR[i][j] -> day i, route j, gives the index in the chromT of the first customer in this route
	vector < vector<int> > chromR ;

	// follows[client][day] gives the next customer in the considered day
	// Used to compute the Hamming distance between solutions
	// if the customer does not exist in this day, returns -1
	vector < vector<int> > follows ;

	// follows[client][day] gives the previous customer in the considered day
	// Used to compute the Hamming distance between solutions
	// if the customer does not exist in this day, returns -1
	vector < vector<int> > previous ;

	// computing the follows and previous tables
	void computeFollows();

	// working table for split (dynamic programming for Split)
	vector < vector < CostSol> > potentials ;

	// table of predecessors (dynamic programming for Split)
	vector < vector < vector<int> > > pred ;

	// storing the cost of the evaluations of arcs (i,j) in the Split graph
	vector < vector < CostSol> > coutArcsSplit ;

	// tells if the fitness has been computed
	bool isFitnessComputed ;

	// tells if the individual is feasible
	bool isValid ;

	// vector used in crossover PIX
	vector < int > toPlace ;

	// measure of distance from "this" to an individual indiv2
	double distance(Individual * indiv2);

	// Individuals ranked by proximity in the population
	list <proxData> plusProches ;

	// functions to manage the "plusProches" structure
	void addProche(Individual * indiv) ;
	void removeProche(Individual * indiv) ;

	// average distance with the n closest individuals
	double distPlusProche(int n) ;

	// Data structure to perform a LocalSearch. 
	// Only some complete individuals, "offspring" for example in Generic.h
	// possess this structure. The others in Population.h are mainly used as containers
	LocalSearch * localSearch ;

	// Data structure for preprocessing information on sequences during the Split algorithm
	vector<SeqData *> seq ;
	SeqData * myseq ;

	// Split function
	// tries first the simple Split without considering the limit on the number of vehicles
	// if the solution does not respect the number of trips, calls the Split with limited fleet.
	void generalSplit();

	// simple Split function based on a shortest path (does not necessary respect the number of vehicles)
	// returns 1 in case of success, 0 otherwise
	int splitSimple(int k) ;

	// Split function for problems with limited fleet, based on a m-shortest path. 
	// (repeating m times the iteration of Bellman on all nodes)
	void splitLF(int k) ;

	// quick function to fill correctly all data structures, once Split has been performed
	void measureSol() ;

	// initialization of the potentials vector for Split
	void initPot(int day) ;

	// updating the LocalSearch structure with the information of the Individual.
	// Warning, Split must have been computed before 
	void updateLS() ;

	// Updating the individual data structures from the local search information
	void updateIndiv() ;

	// little test for debugging (PVRP)
	void testPatternCorrectness();

	// copy of an individual in the other
	// Warning, only copies the chromosomes for storage, not all other structures (potentials and LS)
	void recopieIndividu (Individual * destination , Individual * source);

	// shaking operator, acting on the chromT structure, used by the ILS version of the code
	// with nbShak random swaps between two customer visits in randomly chosen days.
	void shakingSwap (int nbShak);

	// constructor of a random individual
	// if the flag "createAllStructures" is set to true, all search structures, including the LS are also initialized
	Individual(Params * params, bool createAllStructures);

	//destructor
	~Individual();

};
#endif
