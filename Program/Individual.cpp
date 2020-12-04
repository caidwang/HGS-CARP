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

#include "Individual.h"

Individual::Individual(Params * params, bool createAllStructures) : params(params)
{
	// Initializing some structures that enable to compute the data on subsequences in the Split algorithm
	myseq = new SeqData(params);
	for (int i=0 ; i<params->nbClients + params->nbDepots +1 ; i++) 
		seq.push_back(new SeqData(params)); 

	vector <int> tempVect ;
	vector <double> tempVectDbl ;
	pattern p1 ;
	p1.dep = 0 ;
	p1.pat = 0 ;
	localSearch = new LocalSearch() ;

	// Initializing the chromosome structures
	for (int i = 0 ; i <= params->nbDays ; i++)
	{
		chromT.push_back(tempVect);
		chromR.push_back(tempVect);
		for (int j = 0 ; j < params->nombreVehicules[i] ; j++ )
			chromR[i].push_back(-1);
	}

	for (int i = 0 ; i < params->nbClients + params->nbDepots ; i++)
		chromP.push_back(p1);

	for (int i=0 ; i< params->nbDepots + params->nbClients ; i++)
	{
		follows.emplace_back();
		previous.emplace_back();
		for (int k=0 ; k <= params->nbDays ; k++)
		{
			follows[i].push_back(-1);
			previous[i].push_back(-1);
		}
	}

	// If we wish to also create the individual, the local search and Split structures
	if (createAllStructures)
	{
		/* CREATING THE INITIAL INDIVIDUAL */
		vector < vector < int > > tempVect2 ;
		int maxDepots = 0 ;
		int temp, temp2, jj, dayCombinaison, depot  ; 
		age = 0 ;
		isFitnessComputed = false ;
		vector <double> charge ;

		for (int i = 0 ; i <= params->nbDays ; i++)
			charge.push_back(0);

		for (int i = params->nbDepots ; i < params->nbClients + params->nbDepots ; i++)
		{
			if (params->cli[i].freq != 0)
			{
				chromP[i]= params->cli[i].visits[rand() % (int)params->cli[i].visits.size()];
				dayCombinaison = chromP[i].pat ;
				depot = chromP[i].dep ;
				if (chromP[i].dep == -1) cout << "error" << endl ;
				for (int k = 0 ; k < params->formerNbDays ; k++)
				{
					temp2 = dayCombinaison % 2 ;
					dayCombinaison = dayCombinaison/2 ; 
					if (temp2 == 1) chromT[params->formerNbDays - k + depot * params->formerNbDays].push_back(i);
				}
			}
		}

		// shuffling
		for (int k = 1 ; k <= params->nbDays ; k++)
		{
			for (int i = 0 ; i <= (int)chromT[k].size() - 1 ; i++)
			{
				jj = i + rand() % ((int)chromT[k].size() - i) ;
				temp = chromT[k][i] ;
				chromT[k][i] = chromT[k][jj] ;
				chromT[k][jj] = temp ;
			}
		}

		/* CREATING THE SPLIT STRUCTURES */

		for (int k = 0 ; k <= params->nbDays ; k++)
		{
			pred.push_back(tempVect2);
			for (int i = 0 ; i < params->nombreVehicules[k] + 1; i++)
			{
				pred[k].push_back(tempVect);
				pred[k][i].push_back(0);
				for (int j = 0 ; j < (int) params->nbClients + params->nbDepots + 1 ; j++)
					pred[k][i].push_back(0);
			}
		}

		vector <CostSol> potTemp ;
		CostSol csol ;
		csol.evaluation = 1.e30 ;

		for (int k = 0 ; k <= params->nbDays ; k++)
			if (params->nombreVehicules[k] > maxDepots)
				maxDepots = params->nombreVehicules[k] ;

		for (int i = 0 ; i < maxDepots + 1 ; i++)
		{
			potentials.push_back(potTemp);
			for (int j = 0 ; j < (int) params->nbClients + params->nbDepots + 1 ; j++)
				potentials[i].push_back(csol);
		}
        potentials[0][0].evaluation = 0 ;

		for (int i=0 ;  i< params->nbDepots + params->nbClients ; i++)
		{
			coutArcsSplit.push_back(vector<CostSol>());
			for (int j=0 ;  j< params->nbDepots + params->nbClients + 1 ; j++)
				coutArcsSplit[i].push_back(CostSol());
		}
	}
}

void Individual::recopieIndividu (Individual * destination , Individual * source)
{
	destination->chromT = source->chromT ;
	destination->chromP = source->chromP ;
	destination->chromR = source->chromR ;
	destination->costSol.capacityViol = source->costSol.capacityViol ;
	destination->costSol.evaluation = source->costSol.evaluation ;
	destination->costSol.distance = source->costSol.distance ;
	destination->costSol.lengthViol = source->costSol.lengthViol ;
	destination->costSol.routes = source->costSol.routes ;
	destination->isFitnessComputed = source->isFitnessComputed ;
	destination->age = 0 ;
	destination->isValid = source->isValid ;
	destination->follows = source->follows ;
	destination->previous = source->previous ;
	destination->nbRoutes = source->nbRoutes ;
	destination->maxRoute = source->maxRoute ;
	destination->potentials = source->potentials ;
	destination->pred = source->pred ;
	destination->toPlace.clear();
	destination->toPlace = source->toPlace ;
}

void Individual::shakingSwap (int nbShak)
{
	// only used in the ILS
	int itShak = 0 ;
	while (itShak < nbShak)
	{
		// picking a random day to operate the shaking
		int day = rand() % params->nbDays + 1 ;

		// picking two random services in this day to be swapped
		int nbCustDay = (int)chromT[day].size() ;
		int pos1 = rand() % nbCustDay ;
		int pos2 = rand() % nbCustDay ;

		int temp = chromT[day][pos1] ;
		chromT[day][pos1] = chromT[day][pos2] ;
		chromT[day][pos2] = temp ;

		itShak ++ ;
	}
}

Individual::~Individual()
{ 
	FreeClear(seq);
	delete localSearch ;
	delete myseq ;
}

void Individual::generalSplit()
{
    costSol.evaluation = 0 ;
    costSol.capacityViol = 0 ;
    costSol.distance = 0 ;
    costSol.lengthViol = 0 ;
    costSol.routes = 0 ;

	// performing the Split for each day
	// we first try the simple split, 
	for (int k = 1 ; k <= params->nbDays ; k++)
		if (chromT[k].size() != 0  && splitSimple(k) == 0) splitLF(k); // if its not successful we run the Split with limited fleet

	// Do we have a feasible solution
	if (costSol.capacityViol <= 0.00000001 && costSol.lengthViol <= 0.00000001)
        isValid = true ;
	else
        isValid = false;

	// If split was unsuccessful, we allow more capacity violations
	// Usually, the allowed capacity violation is largely sufficient, and this message would indicate a bug.
	if (costSol.evaluation > 1.e20 )
	{
		cout << "Increasing the capacity violation limit in Split" << endl ;
		params->borne *= 1.1 ;
		generalSplit();
	}

	// A quick test for safety (to avoid any possibility of infinite loops and printouts)
	if (params->borne >= 100.0)
		throw string ("Impossible to Split, most likely a problem of the dataset, aborting the run"); 

	isFitnessComputed = true ;
	measureSol(); // calling a post-processing function which fills all other data structures and verifies the solution cost
    computeFollows(); // updating the predecessor and successor structure
}

// Simple Split, does not necessarily respect the number of vehicles
int Individual::splitSimple(int k)
{
	// We will only use the line "0" of the potential and pred data structures
	myseq->initialisation(params->ordreVehicules[k][0].depotNumber,params,this,k,false);
	double cost, mydist,mytminex,myloadex;
	int j ;

	// Main code of Split
	for (int i=0 ; i < (int)chromT[k].size() ; i++ )
	{
		// Compute a route with a single visit
		seq[i]->initialisation(params->ordreVehicules[k][0].depotNumber,params,this,k,false);
		j = i ;
		while (j < (int)chromT[k].size() && seq[j]->load <= params->ordreVehicules[k][0].vehicleCapacity*params->borne )
		{
			// Extend this route to the next visit
			// 这里不会随着i增加被累加， 每次下一个等于前一个seq加上vour，第i个seq永远是0->0，因此可以被重置
			seq[j+1]->concatOneAfter(seq[j],chromT[k][j],this,k);
			cost = seq[0]->evaluation(seq[j+1],myseq,&params->ordreVehicules[k][0],mydist,mytminex,myloadex); // and evaluate

			// Test if this label is better
			if (potentials[0][i].evaluation + cost < potentials[0][j + 1].evaluation )
			{
                potentials[0][j + 1].evaluation = potentials[0][i].evaluation + cost ;
                potentials[0][j + 1].capacityViol = potentials[0][i].capacityViol + myloadex ;
                potentials[0][j + 1].distance = potentials[0][i].distance + mydist ;
                potentials[0][j + 1].lengthViol = potentials[0][i].lengthViol + mytminex ;
                potentials[0][j + 1].routes = potentials[0][i].routes + 1 ;
				pred[k][0][j+1] = i ;
			}
			j++ ;

		}
	}

	// Count the number of routes and see if the soltion is OK
	j = (int)chromT[k].size() ;
	for (int jj = 0 ; jj < params->nombreVehicules[k] ; jj ++ )
	{
		pred[k][params->nombreVehicules[k] - jj][j] = pred[k][0][j] ;
		j = pred[k][params->nombreVehicules[k] - jj][j] ;
	}

	// If we arrived to the beginning
	if (j == 0)
	{
		// We cumulate the costs
		costSol.capacityViol += potentials[0][chromT[k].size()].capacityViol ;
        costSol.evaluation += potentials[0][chromT[k].size()].evaluation ;
        costSol.distance += potentials[0][chromT[k].size()].distance ;
        costSol.lengthViol += potentials[0][chromT[k].size()].lengthViol ;
        costSol.routes += potentials[0][chromT[k].size()].routes ;
		initPot(k); // we reinitialize the dynamic programming structures for the next use
		return 1 ;
	}
	else
	{
		initPot(k); // we reinitialize the dynamic programming structures for the next use
		return 0 ;
	}
}

// Split for problems with limited fleet
// To avoid evaluating several time ("m" times, where m is the number of vehicles) the costs of the arcs 
// (this could be time consuming for some complex VRP variants, for example problems with HOS regulations)
// We already pre-compute and store them in "coutArcsSplit"
void Individual::splitLF(int k)
{ 
	double cost, mydist,mytminex,myloadex;
	int i,j ;
	myseq->initialisation(params->ordreVehicules[k][0].depotNumber,params,this,k,false);

	// preprocessing arc costs
	for (int i=0 ; i < (int)chromT[k].size() ; i++ )
	{
		seq[i]->initialisation(params->ordreVehicules[k][0].depotNumber,params,this,k,false);
		coutArcsSplit[i][i].evaluation = seq[i]->evaluation(seq[i],myseq,&params->ordreVehicules[k][0],mydist,mytminex,myloadex);
		coutArcsSplit[i][i].capacityViol =  myloadex ;
		coutArcsSplit[i][i].distance = mydist ;
		coutArcsSplit[i][i].lengthViol  =  mytminex ;
		for (int j=i ; j < (int)chromT[k].size() && seq[j]->load <= params->ordreVehicules[k][0].vehicleCapacity*params->borne ; j++ )
		{
			seq[j+1]->concatOneAfter(seq[j],chromT[k][j],this,k);
			coutArcsSplit[i][j+1].evaluation = seq[j+1]->evaluation(seq[j+1],myseq,&params->ordreVehicules[k][0],mydist,mytminex,myloadex);
			coutArcsSplit[i][j+1].capacityViol =  myloadex ;
			coutArcsSplit[i][j+1].distance = mydist ;
			coutArcsSplit[i][j+1].lengthViol  =  mytminex ;
		}
	}


	// for each vehicle
	for ( int cam = 0 ; cam < params->nombreVehicules[k] ; cam++)
	{
		i = 0 ;
		// propagate all labels
		while (i < (int)chromT[k].size() && potentials[cam][i].evaluation < 1.e29 )
		{
			cost = coutArcsSplit[i][i].evaluation ;
			myloadex = coutArcsSplit[i][i].capacityViol ;
			mydist =  coutArcsSplit[i][i].distance ;
			mytminex = coutArcsSplit[i][i].lengthViol ;
			if (potentials[cam][i].evaluation + cost < potentials[cam + 1][i].evaluation )
			{
                potentials[cam + 1][i].evaluation = potentials[cam][i].evaluation + cost ;
                potentials[cam + 1][i].capacityViol = potentials[cam][i].capacityViol + myloadex ;
                potentials[cam + 1][i].distance = potentials[cam][i].distance + mydist ;
                potentials[cam + 1][i].lengthViol = potentials[cam][i].lengthViol + mytminex ;
                potentials[cam + 1][i].routes = potentials[cam][i].routes ;
				pred[k][cam+1][i] = i ;
			}
			j = i ;
			while (j < (int)chromT[k].size() && myloadex <= params->ordreVehicules[k][cam].vehicleCapacity*(params->borne-1.0) )
			{
				cost = coutArcsSplit[i][j+1].evaluation ;
				myloadex = coutArcsSplit[i][j+1].capacityViol ;
				mydist =  coutArcsSplit[i][j+1].distance ;
				mytminex = coutArcsSplit[i][j+1].lengthViol ;
				if (potentials[cam][i].evaluation + cost < potentials[cam + 1][j + 1].evaluation )
				{
                    potentials[cam + 1][j + 1].evaluation = potentials[cam][i].evaluation + cost ;
                    potentials[cam + 1][j + 1].capacityViol = potentials[cam][i].capacityViol + myloadex ;
                    potentials[cam + 1][j + 1].distance = potentials[cam][i].distance + mydist ;
                    potentials[cam + 1][j + 1].lengthViol = potentials[cam][i].lengthViol + mytminex ;
                    potentials[cam + 1][j + 1].routes = potentials[cam][i].routes + 1 ;
					pred[k][cam+1][j+1] = i ;
				}
				j++ ;
			}
			i++ ;
		}
	}

	// we cumulate the costs
	costSol.capacityViol += potentials[params->nombreVehicules[k]][chromT[k].size()].capacityViol ;
    costSol.evaluation += potentials[params->nombreVehicules[k]][chromT[k].size()].evaluation ;
    costSol.distance += potentials[params->nombreVehicules[k]][chromT[k].size()].distance ;
    costSol.lengthViol += potentials[params->nombreVehicules[k]][chromT[k].size()].lengthViol ;
    costSol.routes += potentials[params->nombreVehicules[k]][chromT[k].size()].routes ;

	// and clean the dynamic programming structures
	initPot(k);
}

void Individual::measureSol()
{
	int j ;
	nbRoutes = 0 ;
	maxRoute = 0 ;
	double myCost ;
	double totalCost = 0 ;

	for (int kk = 1 ; kk <= params->nbDays ; kk++)
	{
		// we use the result of Split to fill the other data structures
		j = (int)chromT[kk].size() ;
		for (int jj = 0 ; jj < params->nombreVehicules[kk] ; jj ++ )
		{
			if ( j != (int)pred[kk][params->nombreVehicules[kk] - jj][j] ) 
			{
				// Beginning and end of this route in the chromosome
				int deb = (int)pred[kk][params->nombreVehicules[kk] - jj][j] ;
				int end = j-1 ;

				// Counting the number of routes
				nbRoutes++ ; 
				
				// Constructing the data to evaluate the route
				seq[deb]->initialisation(params->ordreVehicules[kk][0].depotNumber,params,this,kk,false);
				for (int i=deb; i <= end ; i++)
					seq[i+1]->concatOneAfter(seq[i],chromT[kk][i],this,kk);
				myCost = seq[deb]->evaluation(seq[end+1],seq[deb],&params->ordreVehicules[kk][0]);
				
				// Measuring this route to see if its longer than the maxRoute
				if (myCost > maxRoute)
					maxRoute = myCost ;

				// And a quick verification of the solution by summing again the cost (for debugging)
				totalCost += myCost ;
			}
			
			j = (int)pred[kk][params->nombreVehicules[kk] - jj][j] ;
			chromR[kk][params->nombreVehicules[kk] - jj - 1] = j ;
		}
	}

	// To avoid any problem of numerical imprecision
	costSol.evaluation = costSol.distance + params->penalityCapa * costSol.capacityViol
                         + params->penalityLength * costSol.lengthViol ;
}

void Individual::initPot(int day)
{
	for (int i = 0 ; i < params->nombreVehicules[day] + 1 ; i++)
		for (size_t j = 0 ; j <= chromT[day].size() + 1 ; j++)
            potentials[i][j].evaluation = 1.e30 ;

    potentials[0][0].evaluation = 0 ;
    potentials[0][0].capacityViol = 0 ;
    potentials[0][0].distance = 0 ;
    potentials[0][0].lengthViol = 0 ;
    potentials[0][0].routes = 0 ;
}

void Individual::updateLS()
{
	// Loading the local search structures
	// Warning, Split must have been computed before
	int depot ;
	int i,j ;
	Node * myDepot ;
	Node * myDepotFin ;
	Node * myClient ;
	Route * myRoute ;

	for (int kk = 1 ; kk <= params->nbDays ; kk++)
	{
        // we reinitialize the "routeOrder" vector
        localSearch->routeOrder[kk].clear();

        // we reset the "isPresent" values to false
        // isPresent 当前客户节点 今天是否被服务
        for (i = params->nbDepots; i < params->nbClients + params->nbDepots; i++)
            localSearch->clients[kk][i].isPresent = false;

        // using Split results to load the solution in the LS structure
        j = (int) chromT[kk].size();

        for (int jj = 0; jj < params->nombreVehicules[kk]; jj++)
		{
			depot = params->ordreVehicules[kk][params->nombreVehicules[kk] - jj - 1].depotNumber ;
			// i -> 前驱节点
			i = (int)pred[kk][params->nombreVehicules[kk] - jj][j] ;

			myDepot = &localSearch->depots[kk][params->nombreVehicules[kk] - jj - 1];
			myDepotFin = &localSearch->depotsFin[kk][params->nombreVehicules[kk] - jj - 1];
			myRoute = &localSearch->routes[kk][params->nombreVehicules[kk] - jj - 1];

			myDepot->nextNode = myDepotFin ;
			myDepot->pred = myDepotFin ;
			myDepotFin->nextNode = myDepot ;
			myDepotFin->pred = myDepot ;

			// a single visit
			if ( j == i+1 ) {
                myClient = &localSearch->clients[kk][chromT[kk][i]];
                myClient->pred = myDepot;
                myClient->nextNode = myDepotFin;
                myClient->route = myRoute;
                myClient->isPresent = true;
                localSearch->nodeTestedForEachRoute(myClient->cour, kk);
                myDepot->nextNode = myClient;
                myDepotFin->pred = myClient;
                localSearch->routeOrder[kk].push_back(myClient->cour);
            }
			else if (j > i+1) {
                // at least two visits
                myClient = &localSearch->clients[kk][chromT[kk][i]];
                myClient->pred = myDepot;
                myClient->nextNode = &localSearch->clients[kk][chromT[kk][i + 1]];
                myClient->route = myRoute;
                myClient->isPresent = true;
                localSearch->nodeTestedForEachRoute(myClient->cour, kk);
                myDepot->nextNode = myClient;
                localSearch->routeOrder[kk].push_back(myClient->cour);

                // information for the end of the route
                myClient = &localSearch->clients[kk][chromT[kk][j - 1]];
                myClient->pred = &localSearch->clients[kk][chromT[kk][j - 2]];
                myClient->nextNode = myDepotFin;
                myClient->route = myRoute;
                myClient->isPresent = true;
                localSearch->nodeTestedForEachRoute(myClient->cour, kk);
                myDepotFin->pred = myClient;
                localSearch->routeOrder[kk].push_back(myClient->cour);

                // and the middle
                for (int k = (int) i + 1; k <= j - 2; k++) {
                    myClient = &localSearch->clients[kk][chromT[kk][k]];
                    myClient->pred = &localSearch->clients[kk][chromT[kk][k - 1]];
                    myClient->nextNode = &localSearch->clients[kk][chromT[kk][k + 1]];
                    myClient->route = myRoute;
                    myClient->isPresent = true;
                    localSearch->nodeTestedForEachRoute(myClient->cour, kk);
                    localSearch->routeOrder[kk].push_back(myClient->cour);
                }
			}
			j = i ;
		}
	
		for (i = 0 ; i < params->nombreVehicules[kk] ; i++ )
			localSearch->routes[kk][i].updateRouteData(false);
	}

	// and we preprocess the route vide structure
	for (int day = 1 ; day <= params->nbDays ; day++)
		localSearch->setRouteVide(day);
}


void Individual::updateIndiv()
{
	// Now, we go through the LS structure to update the individual (its chromosomes)
	int pos ; 
	vector < Route * > ordreRoutes ;
	Node * node ;

	for (int kk = 1 ; kk <= params->nbDays ; kk++)
	{
		ordreRoutes.clear();
		for (int ii=0 ; ii < params->nombreVehicules[kk] ; ii++)
			ordreRoutes.push_back(&localSearch->routes[kk][ii]);

		pos = 0 ;
		for (int r=0 ; r < (int)ordreRoutes.size() ; r ++)
		{
			node = ordreRoutes[r]->depot->nextNode ;
			while (!node->isDepot)
			{
				chromT[kk][pos]= node->cour;
				node = node->nextNode ;
				pos ++ ;
			}
		}
	}

	// And we apply Split to derive all other structures
	// testPatternCorrectness(); // can be used for debugging
	generalSplit();
}

void Individual::testPatternCorrectness()
{
	// Little test for debugging (test that all visits are correct for each customer)
	vector <int> frequencies ;
	for (int i = 0 ; i< params->nbClients + params->nbDepots ; i++)
		frequencies.push_back(0);

	for (int k=1 ; k <= params->nbDays ; k++)
		for (int i=0 ; i<(int)chromT[k].size() ; i++)
			frequencies[chromT[k][i]] ++ ;

	// Here a little test that verifies that all patterns in the chromP correspond to the values in the chromT
	for (int k=1 ; k <= params->nbDays ; k++)
	{
		int realDay = (k-1)%params->formerNbDays + 1 ;
		for (int i=0 ; i<(int)chromT[k].size() ; i++)
		{
			int calcul = chromP[chromT[k][i]].pat ;
			for (int kk = 0 ; kk < params->formerNbDays - realDay ; kk++)
				calcul = calcul/2 ;
			if ( calcul % 2 != 1)
				cout << "Issue, some customers are not placed in accordance to their pattern" << endl ;
		}
	}

	for (int i = params->nbDepots ; i< params->nbClients + params->nbDepots ; i++)
		if (frequencies[i] != params->cli[i].freq) throw string ("Incorrect solution with missing visits") ;

	for (int i = params->nbDepots ; i< params->nbClients + params->nbDepots ; i++)
		if (params->cli[i].visitsDyn[chromP[i].pat] == -1) throw string ("Incorrect solution with missing visits") ;

	for (int i=params->nbDepots ; i < params->nbClients + params->nbDepots ; i++ )
		if ( chromP[i].dep < 0 || chromP[i].dep >= params->nbDepots  ) throw string ("Incorrect solution with missing visits") ;
}

// 计算两个个体之间的hamming距离, 用于判断对多样性的贡献
double Individual::distance(Individual * indiv2)
{
	// Hamming distance
	bool isIdentical ;
	double note = 0.0 ;

	for (int j=params->nbDepots ; j < params->nbClients + params->nbDepots ; j++)
	{
		// For PVRP (Hamming distance based on the patterns)
		// 和文章描述有出入, 文章中, 如果pattern不同并且depot不同, 应该是累加2, 此处累加1
		isIdentical = (chromP[j].pat == indiv2->chromP[j].pat && chromP[j].dep == indiv2->chromP[j].dep) ;

		// For CVRP (Hamming distance based on the predecessors/successors)
		if (!params->periodique)
		{
			for (int s=1  ; s <= params->nbDays ; s++)
			{
				if ((follows[j][s] != indiv2->follows[j][s] || previous[j][s] != indiv2->previous[j][s] )
					&& (previous[j][s] != indiv2->follows[j][s] || follows[j][s] != indiv2->previous[j][s] ))
					isIdentical = false ;
			}
		}

		if (!isIdentical)
			note += 1.0 ;
	}

	return ((double)note /(double)(2*params->nbClients)) ;
}

// updating the predecessor and successor structure
void Individual::computeFollows ()
{
	int jj ;
	for (int i=0 ; i< params->nbDepots + params->nbClients ; i++)
	{
		for (int k=1 ; k<= params->nbDays; k++)
		{
            follows[i][k] = -1 ;
            previous[i][k] = -1 ;
		}
	}

	for (int k=1 ; k <= params->nbDays ; k++)
	{
		if (chromT[k].size() != 0)
		{
			for (int i=0 ; i < (int)chromT[k].size()-1 ; i++)
                follows[chromT[k][i]][k] = chromT[k][i + 1];

			for (int i=1 ; i < (int)chromT[k].size() ; i++)
                previous[chromT[k][i]][k] = chromT[k][i - 1];

            follows[chromT[k][chromT[k].size() - 1]][k] = params->ordreVehicules[k][0].depotNumber;
            previous[chromT[k][0]][k] = params->ordreVehicules[k][0].depotNumber;

			// arranging those which are located at the beginning or end of a route
			for (int i=0 ; i < params->nombreVehicules[k] ; i++)
			{
				jj = chromR[k][i] ;
                previous[chromT[k][jj]][k] = params->ordreVehicules[k][0].depotNumber;
				if (jj != 0)
                    follows[chromT[k][jj - 1]][k] = params->ordreVehicules[k][0].depotNumber;
			}
		}
	}
}

// 维护个体在种群中相似个体的列表, 列表中个体与该个体的相似程度逐渐降低 proche -> close 相似
void Individual::addClose(Individual *indiv) {
    // Adding an individual in the structure of proximity (diversity management procedures)
    list<proxData>::iterator it;
    proxData data; // data about the proximity of an Individual with regards to the others in the population.
    data.indiv = indiv;
    data.dist = distance(indiv);

    if (plusProches.empty()) plusProches.push_back(data);
    else {
		it = plusProches.begin();
		while ( it != plusProches.end() && it->dist < data.dist)
			++it ;
		plusProches.insert(it,data);
	}
}

void Individual::removeProche(Individual * indiv)
{
	// Removing an individual in the structure of proximity (diversity management procedures)
	list<proxData>::iterator last = plusProches.end();
	for (list<proxData>::iterator first = plusProches.begin() ; first != last ; )
		if (first->indiv == indiv)
			first = plusProches.erase(first) ;
		else
			++first ;
}

double Individual::distPlusProche(int n)
{
	// Computing the average distance with the close elements (diversity management)
	double result = 0 ;
	double compte = 0 ;
	list<proxData>::iterator it = plusProches.begin();

	for (int i=0 ; i<n && it!= plusProches.end(); i++)
	{
		result += it->dist ;
		compte += 1.0 ;
		++it ;
	}
	return result/compte ;
}
