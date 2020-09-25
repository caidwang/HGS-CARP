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

#include "Genetic.h"

void Genetic::evolve (int maxIterNonProd, int nbRec)
{
	// Should we run either ILS or and HGA 
	if (params->isILS_general) 
		evolveILS();
	else 
		evolveHGA(maxIterNonProd, nbRec);
}

void Genetic::evolveHGA (int maxIterNonProd, int nbRec)
{
	// Main code of the HGA 
	Individual * parent1 ;
	Individual * parent2 ;
	int place, place2 ;
	place2 = 10000 ;
	nbIterNonProd = 1 ;
	nbIter = 0 ;
	int resultCross ;
	string temp ;
	double fitBeforeRepair ;
	CostSol bestSolFeasibility ;
	clock_t start = clock() ; // When iterating several time the HGA (e.g. PCARP, the time limit applies to one iteration -- fleet size or max distance value)

	if (population->getIndividuBestValide() != NULL) bestSolFeasibility = population->getIndividuBestValide()->costSol ;
	else bestSolFeasibility = population->getIndividuBestInvalide()->costSol ;
	for (int i=0 ; i<population->invalides->nbIndiv ; i++)
		if (population->invalides->individus[i]->costSol.isBetterFeas(bestSolFeasibility)) bestSolFeasibility = population->invalides->individus[i]->costSol ;

    offspring->localSearch->nbTotalRISinceBeginning = 0 ;
    offspring->localSearch->nbTotalPISinceBeginning = 0 ;

	cout << "| Start of GA | NbNodes : " << params->nbClients << " | NbVehicles : " << params->nbVehiculesPerDep << " | " << endl ;

	while (nbIterNonProd < maxIterNonProd && (clock() - start <= ticks) && (!params->isSearchingFeasible || population->getIndividuBestValide() == NULL))
	{
		// CROSSOVER
		parent1 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		parent2 = population->getIndividuBinT(); // Pick two individuals per binary tournament
		offspring->recopieIndividu(offspring, parent1); // Put them in adequate data structures
		offspring2->recopieIndividu(offspring2, parent2); // Put them in adequate data structures

		if (!params->periodique && !params->multiDepot) 
			resultCross = crossOX(); // Pick OX crossover if its a single-period problem
		else 
			resultCross = crossPIX() ; // Otherwise PIX (see Vidal et al 2012 -- OR)

		// SPLIT
		offspring->generalSplit();

		// LOCAL SEARCH
		offspring->updateLS();
		offspring->localSearch->runSearchTotal();
		offspring->updateIndiv();
		population->updateNbValides(offspring);
		place = population->addIndividu(offspring) ;

		// POSSIBLE REPAIR
		if (!offspring->isValid)
		{
			fitBeforeRepair = offspring->costSol.evaluation ;
			if (rand() % 2 == 0) // 50% chance to do repair on an infeasible individual
			{
                repair();
				if (offspring->costSol.evaluation < fitBeforeRepair - 0.01 || offspring->costSol.evaluation > fitBeforeRepair + 0.01 || offspring->isValid)
					place2 = population->addIndividu(offspring) ;
				if (offspring->isValid)
					place = place2 ;
				else 
					place = min(place,place2);
			}
		}

		// SOME TRACES
		if ((offspring->isValid && place == 0) || (offspring->costSol.isBetterFeas(bestSolFeasibility) && population->valides->nbIndiv == 0))
		{	
			if (traces && population->valides->nbIndiv > 0) 
				cout << "NEW BEST FEASIBLE " << place << " " << population->getIndividuBestValide()->costSol.evaluation << " distance : " << offspring->costSol.distance << " nbRoutes : " << offspring->costSol.routes << " capaViol : " << offspring->costSol.capacityViol << " lengthViol : " << offspring->costSol.lengthViol << endl << endl ;
			if (traces && population->valides->nbIndiv == 0 ) 
				cout << "NEW BEST INFEASIBLE " << place << " " << offspring->costSol.evaluation << " distance : " << offspring->costSol.distance << " nbRoutes : " << offspring->costSol.routes << " capaViol : " << offspring->costSol.capacityViol << " lengthViol : " << offspring->costSol.lengthViol << endl << endl ;
			if (offspring->costSol.isBetterFeas(bestSolFeasibility))
				bestSolFeasibility = offspring->costSol ;
			nbIterNonProd = 1 ; 
		}
		else nbIterNonProd ++ ;


		// DIVERSIFICATION
		if (nbRec > 0 && nbIterNonProd % (maxIterNonProd/3+1) == maxIterNonProd/3) 
		{	
			if (traces) cout << "Diversification" << endl ;
			population->diversify();
		}

		// PENALTY MANAGEMENT
		if (nbIter % 30 == 0) 
			gererPenalites () ;

		// MORE TRACES
		if (traces && nbIter % 500 == 0)
		{
			population->afficheEtat(nbIter);
			cout << " | NbTotalMovesLS : (RI) " << offspring->localSearch->nbTotalRISinceBeginning << " | (PI) " << offspring->localSearch->nbTotalPISinceBeginning << " | " << endl ;
			cout << " | interSwap " << offspring->localSearch->nbInterSwap ;
			cout << " | intraSwap " << offspring->localSearch->nbIntraSwap ;
			cout << " | inter2opt " << offspring->localSearch->nbInter2Opt ;
			cout << " | intra2opt " << offspring->localSearch->nbIntra2Opt ;
			cout << " | " << endl ;
			cout << " | CPU Time : " <<  (double)clock()/(double)CLOCKS_PER_SEC << " seconds " << endl ;
			cout << endl ;
		}
		nbIter ++ ;
	}

	// END OF THE ALGORITHM
	if (traces)
	{
		cout << "Time Elapsed : " << clock() << endl ;
		cout << "Number of Iterations : " << nbIter << endl ;
	}
}

void Genetic::evolveILS ()
{
	int nbGRASP = 5 ;
	int nbILS = 100 ;
	int nbCHILD = 50 ;
	bool isFirstLoop ;
	Individual * parent ;
	clock_t timeBest2 ;
    offspring->localSearch->nbTotalRISinceBeginning = 0 ;
    offspring->localSearch->nbTotalPISinceBeginning = 0 ;
	nbIter = 0 ;
	clock_t debut = clock();
    offspringBestFoundAll->costSol.evaluation = 1.e30 ;

	cout << "| Debut evolution ILS | NbNodes : " << params->nbClients << " | NbVehicles : " << params->nbVehiculesPerDep << " | " << endl ;

	for (int phaseGrasp = 0 ; phaseGrasp < nbGRASP ; phaseGrasp ++)
	{
		// NEW RANDOM START
		cout << endl  << "------------ RESTART ---------" << endl << endl ;
		params->penalityCapa = 50 ;
		params->penalityLength = 50;
		Individual * myIndiv = new Individual(params, 1.0) ;
		offspringP1->recopieIndividu(offspringP1, myIndiv);
		delete myIndiv ;
        offspringBestFound->costSol.evaluation = 1.e30 ;
		isFirstLoop = true ;
		
		for (int nbGenerationNonProd = 0 ; nbGenerationNonProd < nbILS && (clock() - debut <= ticks) ; nbGenerationNonProd ++)
		{
			if (!isFirstLoop)
			{
				parent = population->getIndividuBestValide();
				if (parent == NULL) parent = population->getIndividuBestInvalide();
				offspringP1->recopieIndividu(offspringP1, parent);
			}
			else 
				isFirstLoop = false ;

			// Clear the population
			population->clear();

			for (int i=0 ; i < nbCHILD ; i++)
			{
				// SHAKING
				offspring->recopieIndividu(offspring, offspringP1);
				offspring->shakingSwap(2 + (int)params->nbClients / 200);
				offspring->generalSplit();

				// LOCAL SEARCH
				offspring->updateLS();
				offspring->localSearch->runSearchTotal();
				offspring->updateIndiv();
				population->updateNbValides(offspring);
				// If the solution is infeasible, do a LS with higher penalty to restore feasibiliy
				if (!offspring->isValid) repair();
				population->addIndividu(offspring) ;

				// Checking if its a solution improvement
				if (offspring->isValid && offspring->costSol.evaluation < offspringBestFound->costSol.evaluation - 0.001)
				{	
					nbGenerationNonProd = -1 ;
					offspringBestFound->recopieIndividu(offspringBestFound, offspring);
					if (offspringBestFound->costSol.evaluation < offspringBestFoundAll->costSol.evaluation - 0.001)
					{
						cout << "NEW BEST EVER : " << offspringBestFound->costSol.evaluation << endl ;
						offspringBestFoundAll->recopieIndividu(offspringBestFoundAll, offspringBestFound);
						timeBest2 = population->timeBest ;
					}
					else
						cout << "NEW BEST      : " << offspring->costSol.evaluation << endl ;
				}

				// Regular adaptation of penalty parameters (every 30 LS)
				if (nbIter % 30 == 0) 
					gererPenalites () ;
				nbIter ++ ;
			}

			if (nbIter % 500 == 0)
			{
				population->afficheEtat(nbIter);
				cout << " | NbTotalMovesLS : (RI) " << offspring->localSearch->nbTotalRISinceBeginning << " | (PI) " << offspring->localSearch->nbTotalPISinceBeginning << " | " << endl ;
				cout << " | interSwap " << offspring->localSearch->nbInterSwap ;
				cout << " | intraSwap " << offspring->localSearch->nbIntraSwap ;
				cout << " | inter2opt " << offspring->localSearch->nbInter2Opt ;
				cout << " | intra2opt " << offspring->localSearch->nbIntra2Opt ;
				cout << " | " << endl ;
				cout << " | CPU Time : " <<  (double)clock()/(double)CLOCKS_PER_SEC << " seconds " << endl ;
				cout << endl ;
			}
		}
	}

	// end of the algorithm, various information displayed
	if (traces) 
	{
		cout << "temps passe : " << clock() << endl ;
		cout << "fin evolution ILS, nombre d'iterations : " << nbIter << endl ;
	}

	// add the best solution found in the population for the writing at the end of the algorithm
	population->addIndividu(offspringBestFoundAll);
	population->timeBest = timeBest2 ; // need to correct the time to best solution.
}

void Genetic::repair ()
{
	double temp, temp2  ;

	temp = params->penalityCapa ;
	temp2 = params->penalityLength ;

	// First Tentative of Repair
	params->penalityCapa *= 10 ;
	params->penalityLength *= 10 ;
	offspring->updateLS();
	offspring->localSearch->runSearchTotal();
	offspring->updateIndiv();

	// If the first tentative failed, second tentative with higher penalty
	if (!offspring->isValid)
	{
		params->penalityCapa *= 10 ;
		params->penalityLength *= 10 ;
		offspring->generalSplit();
		offspring->updateLS();
		offspring->localSearch->runSearchTotal();
		offspring->updateIndiv();
	}

	params->penalityCapa = temp ;
	params->penalityLength = temp2 ;
	offspring->measureSol();
}

void Genetic::gererPenalites ()
{
	bool changeDone = false ;
	double fractionCharge = population->fractionValidesCharge() ;
	double fractionTemps = population->fractionValidesTemps() ;

	// if there are not enough feasible solutions
	if ( fractionCharge < params->minValides && params->penalityCapa < 5000)
	{
		params->penalityCapa = (double)((float)params->penalityCapa * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionCharge > params->maxValides && params->penalityCapa > 0.01)
	{
		params->penalityCapa =  (double)((float)params->penalityCapa * 0.85) ;
		changeDone = true ;
	}

	// if there are too many feasible solutions
	if ( fractionTemps < params->minValides && params->penalityLength < 5000)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 1.2) ;
		changeDone = true ;
	}
	else if ( fractionTemps > params->maxValides && params->penalityLength > 0.01)
	{
		params->penalityLength =  (double)((float)params->penalityLength * 0.85) ;
		changeDone = true ;
	}

	if (changeDone) 
		population->validatePen(population->invalides);
}

int Genetic::crossOX ()
{
	int temp, tempSuiv ;

	// We pick the beginning and end of the crossover zone
	int debut = rand() % params->nbClients ;
	int fin = rand() % params->nbClients ;
	while (fin == debut && params->nbClients > 1)
		fin = rand() % params->nbClients ;

	// We initialize a little frequency table to know if each customer was placed or not
	for (int i=params->nbDepots ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient[i] = 1 ;

	int j = debut ;
	// we keep the elements from "debut" to "end"
	while ((j % params->nbClients) != ((fin + 1) % params->nbClients))
	{ 
		freqClient[offspring->chromT[1][j % params->nbClients]] = 0 ;
		j ++ ;
	}

	// We fill the rest of the elements in the order of the second parent
	for (int i=1 ; i <= params->nbClients ; i++)
	{
		temp = offspring2->chromT[1][(fin + i) % params->nbClients] ;
		tempSuiv = offspring2->chromT[1][(fin + i + 1) % params->nbClients] ;
		if (freqClient[temp] == 1)
		{
            offspring->chromT[1][j % params->nbClients] = temp ;
			j ++ ;
		}
	}

	return 0 ;
}

int Genetic::crossPIX ()
{
	vector < int > empty, keep, daysDisturb, boardEnd, boardStatus ; // vide -- 空闲的  garder -- 保存 daysDisturb -- 扰动日期
	vector < vector <int> > keep2 ;
	int jj,i,ii,j,temp,size ;
	int start, end, day ;
	int j1,j2 ;
	int placeInsertion ;

	// We initialize some little data structures
    initialInCrossPIX();

    // We take the days in random order
    disturbDays(daysDisturb);

    // We pick j1 and j2. j1 < j2 < nbDays
    generateN1AndN2ForSplitPopulation(j1, j2);

    // init keep2, store the incomplete chromT
	for (int k = 0 ; k <= params->nbDays ; k ++ )
		keep2.push_back(empty);

	// Inheriting the visits
	for (int k = 0 ; k < params->nbDays ; k ++ )
	{
		day = daysDisturb[k];

		// First case, we copy a segment (these visits will be temporarily kept in the data structure "keep2") ^mix
		if (k < j1 && !offspring->chromT[day].empty())
		{
		    // generate alpha_i and beta_i
            start = (int)(rand() % offspring->chromT[day].size()) ;
            end = (int)(rand() % offspring->chromT[day].size()) ;
			boardEnd.push_back(end);
			j = start ;
			while ( j != ((end + 1) % offspring->chromT[day].size()) )
			{
                updateDepotAndPattern(offspring->chromT[day][j], day);
				keep2[day].push_back(offspring->chromT[day][j]) ;
				j = (j+1) % offspring->chromT[day].size() ;
			}
			offspring->chromT[day].clear();
		}
		else if (k < j2) // Second case, we copy nothing on this day ^2
		{
			offspring->chromT[day].clear() ;
			boardEnd.push_back(-1);
		}
		else // Third case, we copy everything on this day ^1 set
		{
			boardEnd.push_back(0);
			// 遍历的是某一天的所有tour
			for (j=0 ; j < (int)offspring->chromT[day].size() ; j++)
			{
                updateDepotAndPattern(offspring->chromT[day][j], day);
			}
		}
	}

	// We complete with the second parent
	for (int k = 0 ; k < params->nbDays ; k ++ )
	{
		day = daysDisturb[k] ;
        end = boardEnd[k] ;
		if (k < j2)
		{
			for (i=0 ; i < (int)offspring2->chromT[day].size() ; i++)
			{
				ii = offspring2->chromT[day][(i + end + 1) % (int)offspring2->chromT[day].size() ] ;
				if (freqClient[ii] != 0
					&& params->cli[ii].jourSuiv[offspring->chromP[ii].pat][(day - 1) % params->formerNbDays + 1] == (int)((day - 1) % params->formerNbDays + 1)
					&& (offspring->chromP[ii].dep == -1 || offspring->chromP[ii].dep == (day - 1) / params->formerNbDays ))
				{
					offspring->chromT[day].push_back(ii);
                    updateDepotAndPattern(ii, day);
				}
			}
		}
	}

	// we complete with the elements of keep2 (which come from the first parent for the days where the parents are mixed)
	for (int k=1 ; k<=params->nbDays ; k++)
	{
		keep.clear();
		// Choose a random place of insertion
		size = (int)offspring->chromT[k].size() ;
		if (size != 0) placeInsertion = rand() % size ; 
		else placeInsertion = 0 ;
		for (int iii=placeInsertion ; iii <  size ; iii ++)
			keep.push_back(offspring->chromT[k][iii]);
		for (int iii=placeInsertion ; iii <  size ; iii ++)
			offspring->chromT[k].pop_back();
		for (int iii=0 ; iii < (int)keep2[k].size() ; iii ++)
			offspring->chromT[k].push_back(keep2[k][iii]);
		for (int iii=0 ; iii < (int)keep.size() ; iii ++)
			offspring->chromT[k].push_back(keep[iii]);
	}

	offspring->toPlace.clear();
	// We gather in "toPlace" those elements with missing visits
	for (i=0 ; i < params->nbClients + params->nbDepots ; i++ )
	{
		if (freqClient[i] > 0)
		{
			for (j=0 ; j< freqClient[i] ; j++)
				offspring->toPlace.push_back(i);
		}
	}

	// We randomize toPlace
	for (i = 0 ; i < (int)offspring->toPlace.size() ; i++)
	{
		jj = i + rand() % ((int)offspring->toPlace.size() - i) ;
		temp = offspring->toPlace[i] ;
        offspring->toPlace[i] = offspring->toPlace[jj] ;
        offspring->toPlace[jj] = temp ;
	}

	// We call Split to obtain a full solution
	offspring->generalSplit();
	offspring->updateLS();
	offspring->localSearch->placeManquants(); // and we perform a least-cost insertion of the missing visits
	offspring->updateIndiv();

	return 0 ;
}

void Genetic::generateN1AndN2ForSplitPopulation(int &j1, int &j2) const {
    int temp;
    j1 = rand() % params->nbDays ;
    j2 = rand() % params->nbDays ;
    if (j1 > j2)
    {
        temp = j2 ;
        j2 = j1 ;
        j1 = temp ;
    }
}

void Genetic::initialInCrossPIX()  {
    int i;
    for (i=0 ; i < params->nbClients + params->nbDepots ; i++ )
    {
        freqClient[i] = params->cli[i].freq ;
        offspring->chromP[i].pat = 0 ;
        offspring->chromP[i].dep = -1 ;
    }
}

void Genetic::disturbDays(vector<int> &daysDisturb) const {
    int i, jj, temp;
    for (int k=1 ; k <= params->nbDays ; k++ )
        daysDisturb.push_back(k) ;
    for (i = 0 ; i < (int)daysDisturb.size() ; i++)
    {
        jj = i + rand() % ((int)daysDisturb.size() - i) ;
        temp = daysDisturb[i] ;
        daysDisturb[i] = daysDisturb[jj] ;
        daysDisturb[jj] = temp ;
    }
}

void Genetic::updateDepotAndPattern(int customerVertexIdx, int day) {
    freqClient[customerVertexIdx] -= 1 ;
    // pattern 实际上是按位编码的, day0 在最高位 dayN在最低位
    offspring->chromP[customerVertexIdx].pat += (int)pow((float)2, (int)((params->nbDays - day) % params->formerNbDays)) ;
    // dep的赋值: nbDay 实际值是 nbDay * nbDepot 通过 day - 1 / nbDays 落到实际的depot的编号上
    offspring->chromP[customerVertexIdx].dep = (day - 1) / params->formerNbDays ;
}

Genetic::Genetic(Params * params,Population * population, clock_t ticks, bool traces) :
params(params) , population(population) , ticks(ticks) , traces(traces)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i++ )
		freqClient.push_back(params->cli[i].freq);

	// Creating the Individuals that serve to perform the Local Search and other operations
	offspring = new Individual (params, true) ;
    offspring2 = new Individual(params, true) ;
    offspringP1 = new Individual(params, true) ;
    offspringP2 = new Individual(params, true) ;
    offspringBestFound = new Individual(params, true) ;
    offspringBestFoundAll = new Individual(params, true) ;
    offspring->localSearch = new LocalSearch(params, offspring) ;
}

Genetic::~Genetic(void)
{ 
	delete offspring ;
	delete offspring2 ;
	delete offspringP1 ;
	delete offspringP2 ;
	delete offspringBestFound ;
	delete offspringBestFoundAll ;
}

