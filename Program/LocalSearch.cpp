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

#include "LocalSearch.h"
#include "Individual.h"

void LocalSearch::runSearchTotal ()
{
    // Shuffling the order of move evaluations
    params->shuffleClose();
    shuffleRouteOrder();
    int nbMoves = 0;

    /* RUNNING THE LS */
    // RI -- Route improvement
    updateMoves(); // 总体更新的是clients[day][clientIndex]中的moves结构，在判断closestNeighbor时，将该client视为close的，并且在当天需要服务的客户节点，作为待选择的动作
    reinitAllSingleDayMoves();
    for (int day = 1; day <= params->nbDays; day++)
        nbMoves += mutationSameDay(day);
    nbTotalRISinceBeginning += nbMoves;

	// PI and RI local search (only for PCARP and MDCARP, see Vidal et al 2012 (OR)
	if (params->periodique || params->multiDepot)
	{
		// PI -- Pattern improvement
		nbMoves = mutationDifferentDay () ; 
		nbTotalPISinceBeginning += nbMoves ;

		if (nbMoves > 0)
		{
			// RI -- Route improvement
			updateMoves (); 
			reinitAllSingleDayMoves();
			nbMoves = 0 ;
			for (int day = 1 ; day <= params->nbDays ; day++)
				nbMoves += mutationSameDay (day) ;
			nbTotalRISinceBeginning += nbMoves ;
		}
	}
}

// main function for the RI procedure 突变过程
int LocalSearch::mutationSameDay (int day) {
    // Local Search for one given day
    // Based on the classic moves (Relocate, Swap, CROSS, 2-Opt and 2-Opt*)
    researchCompleted = false;
    double costBeforeRemoval, costAfterRemoval;
    bool gainWhenRemoving;
    int movePerformed = 0;
    int nbMoves = 0;
    bool emptyRouteTested;
    Node *tempNode;
    int size2;

    /*
     * 1. 根据routeOrder的顺序，选取U节点
     *    - 如果U节点的移除不能获得收益，之后不惊醒swap和relocate相关的邻域操作
     * 2. 在U节点的所有moves节点中顺序选择节点：
     *    - 如果移除可以获得收益，Relocate, Swap, CROSS and I-CROSS
     *    - 如果V和U不在同一个路径，2-opt 2-opt*
     *    - 将V设为U后的节点，进行路径内的2opt
     * 3. 测试将U插在路径头部和插入一个空路径
     * 4. EjectChain阶段
     */

    // We search and apply moves until a local minimum is attained
    while (!researchCompleted) {
        researchCompleted = true;
        movePerformed = 0;

        // For every node U in random order
        for (int posU = 0; posU < (int) routeOrder[day].size(); posU++) {
            // movePerformed 只能是0/1，当采纳了一次邻域操作时，说明这个节点的变换可能还没到达局部最优，下一轮邻域尝试依然在这个节点上进行
            posU -= movePerformed; // We return on the last node if there has been a move
            nbMoves += movePerformed;
            movePerformed = 0;
            emptyRouteTested = false;

            nodeU = &clients[day][routeOrder[day][posU]];
            routeU = nodeU->route;
            vehicleU = routeU->vehicle;
            x = nodeU->nextNode;

            // In the CARP, some service removal do not reduce the cost of the route
            // In this case, very few moves can improve the solution (only 2-opt variants), and thus SWAP and RELOCATE variants do not need to be tested
            // 去除不必要的LS步骤以缩小LS的搜索空间
            costBeforeRemoval = nodeU->seq0_i->evaluation(routeU->depot->pred->seq0_i, routeU->vehicle); // 去最后一个depot取0->i的seq的penalty cost
            costAfterRemoval = nodeU->seq0_i->evaluation(nodeU->pred->seq0_i, nodeU->nextNode->seqi_n, routeU->vehicle);
            gainWhenRemoving = (costAfterRemoval < costBeforeRemoval - 0.1);

            // For every node V in random order
            size2 = (int) nodeU->moves.size();
            // U,V两个node的操作 O(N^2)
            for (int posV = 0; posV < size2 && movePerformed == 0; posV++) {
                nodeV = &clients[day][nodeU->moves[posV]];
                routeV = nodeV->route;
                vehicleV = routeV->vehicle;

                // If we have not yet tested the moves involving the node U and the route of node V
                // (This flag is reset to false as soon as there is a modification in the route)
                if (!nodeV->route->nodeAndRouteTested[nodeU->cour]) {
                    y = nodeV->nextNode;
                    if (routeV->cour != routeU->cour) {
                        if (movePerformed != 1) {
                            tempNode = nodeV;
                            nodeV = nodeV->nextNode;
                            y = nodeV->nextNode;
                            // Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV
                            // Case where they are in different routes
                            if (gainWhenRemoving) movePerformed = interRouteGeneralInsert();
                            nodeV = tempNode;
                            y = nodeV->nextNode;
                        }

                        // 2-Opt*
                        if (movePerformed != 1)
                            movePerformed = interRoute2Opt();

                        // 2-Opt* (second type, where the routes can be reversed)
                        if (movePerformed != 1)
                            movePerformed = interRoute2OptInv();
                    } else {
                        tempNode = nodeV;
                        nodeV = nodeV->nextNode;
                        y = nodeV->nextNode;

                        // Testing Relocate, Swap, CROSS and I-CROSS (limited to 2 customers) of nodeU and nodeV
                        // Case where they are in the same route
                        if (movePerformed != 1 && gainWhenRemoving)
                            movePerformed = intraRouteGeneralInsertDroite();
                        nodeV = tempNode;
                        y = nodeV->nextNode;
                    }
                }
            }

            nodeV = nodeU->nextNode ;
            routeV = nodeV->route ;
            vehicleV = routeV->vehicle ;
            y = nodeV->nextNode ;
            if (!nodeV->route->nodeAndRouteTested[nodeU->cour]) {
                while (movePerformed != 1 && !nodeV->isDepot) {
                    // Testing 2-Opt between U and V (if the restriction of the granular search allows)
                    if (params->isCorrelated[nodeU->pred->cour][nodeV->cour] ||
                        params->isCorrelated[nodeU->cour][nodeV->nextNode->cour])
                        movePerformed = intraRoute2Opt();
                    nodeV = nodeV->nextNode;
                    y = nodeV->nextNode;
                }
            }

            // Special cases : testing the insertions behind the depot, and the empty routes
            for (int route = 0; route < params->nbVehiclesPerDepot && movePerformed == 0; route++) {
                nodeV = &depots[day][route];
                routeV = nodeV->route;
                y = nodeV->nextNode;
                if ((!nodeV->route->nodeAndRouteTested[nodeU->cour]) && (!y->isDepot || !emptyRouteTested)) {
                    if (y->isDepot) emptyRouteTested = true;
                    if (routeV != routeU) {
                        tempNode = nodeV;
                        nodeV = depots[day][route].nextNode;
                        y = nodeV->nextNode;

                        // Insertion after the depot, in a different route
                        if (gainWhenRemoving && (params->isCorrelated[nodeU->cour][nodeV->cour] ||
                                                 params->isCorrelated[nodeU->cour][y->cour]) && movePerformed != 1)
                            movePerformed = interRouteGeneralInsert();

                        nodeV = nodeV->route->depot;
                        y = nodeV->nextNode;

                        // 2-Opt* after the depot
                        if (params->isCorrelated[nodeU->pred->cour][nodeV->cour] && movePerformed != 1)
                            movePerformed = interRoute2Opt();

                        // 2-Opt* after the depot
                        if ((params->isCorrelated[x->cour][y->cour] || params->isCorrelated[y->cour][x->cour]) &&
                            movePerformed != 1)
                            movePerformed = interRoute2OptInv();

                        nodeV = tempNode;
                        y = nodeV->nextNode;
                    } else {
                        tempNode = nodeV;
                        nodeV = depots[day][route].nextNode;
                        y = nodeV->nextNode;

                        // Insertion after the depot, in the same route
                        if ((params->isCorrelated[nodeU->cour][nodeV->cour] ||
                             params->isCorrelated[nodeU->cour][y->cour]) && movePerformed != 1)
                            movePerformed = intraRouteGeneralInsertDroite();

                        nodeV = tempNode;
                        y = nodeV->nextNode;
                    }
                }
            }

            // Say that we have tested the node U with all routes
            if (movePerformed == 0)
                nodeTestedForEachRoute(nodeU->cour, day);
        }
    }
    // Calling the ejection chains at the end of the LS
    nbMoves += ejectionChains(day);
    return nbMoves ;
}

int LocalSearch::mutationDifferentDay ()
{
	// Local Search to improve the pattern choices for customers (PCARP and MDCARP)
	// Only a single move, which is a relocate of all occurences of a customer in the best combination of days
	// This move is still used for problems with a single period, in this case it only does a relocate without granular search restriction
	researchCompleted = false ;
	int nbMoves = 0 ;
	firstLoop = true ;

	for (int day = 1 ; day <= params->nbDays ; day++)
		for (int r=0 ; r < params->numberVehicle[day] ; r++)
			routes[day][r].initiateInsertions() ;

	while ( !researchCompleted )
	{
        researchCompleted = true ;
		// Searching for a better insertion place for all customers
        for (int posU = 0; posU < params->nbClients; posU++)
            nbMoves += searchBetterPattern(routeOrder[0][posU]);
		firstLoop = false ;
	}
	return nbMoves ;
}

int LocalSearch::interRouteGeneralInsert()
{
	// For a pair of nodes U, V, tests together the Relocate, Swap, and variants of CROSS and I-CROSS limited to two consecutive nodes.
	// Some route evaluations can be gained (about 40%) by doing these computations is a combined manner
	int ibest, jbest ;
	double moveMin ;
	double temp ;
	SeqData * seq = nodeU->seq0_i ;
	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ;

	// This table will keep the result of the moves
	// 0 -> send nothing
	// 1 -> send U
	// 2 -> send U,Unext
	// 3 -> send Unext,U
	for (int i=0 ; i<4 ; i++ )
		for (int j=0 ; j<4 ; j++)
			resultMoves[i][j] = 1.e29 ; 

	// This keeps the current cost
	resultMoves[0][0] = costZero ;

	/* EVALUATION OF MOVE LOWER BOUNDS */
	// We start to compute lower bounds on move evaluations, so that we can discard any moves which are not promising.

	// Can send something if U is not a depot
	if (!nodeU->isDepot)
	{
		// We evaluate the result of the move as a combination of existing subsequences
		resultMoves[1][0] = seq->evaluationLB(nodeU->pred->seq0_i, x->seqi_n, routeU->vehicle);

		// Can we receive something (if V is not a depot)
		// In the following, all the conditionals are set up to avoid moving depots
		if (!nodeV->isDepot)
			resultMoves[1][1] = seq->evaluationLB(nodeU->pred->seq0_i, nodeV->seq1, x->seqi_n, routeU->vehicle);

		// And so on...
		if (!x->isDepot)
		{
			resultMoves[2][0] = seq->evaluationLB(nodeU->pred->seq0_i, x->nextNode->seqi_n, routeU->vehicle);
			resultMoves[3][0] = resultMoves[2][0] ;
			if (!nodeV->isDepot)
			{
				resultMoves[2][1] = seq->evaluationLB(nodeU->pred->seq0_i, nodeV->seq1, x->nextNode->seqi_n, routeU->vehicle);
				resultMoves[3][1] = resultMoves[2][1] ;
				if (!y->isDepot)
				{
					resultMoves[2][2] = seq->evaluationLB(nodeU->pred->seq0_i, nodeV->seq12, x->nextNode->seqi_n, routeU->vehicle);
					resultMoves[2][3] = seq->evaluationLB(nodeU->pred->seq0_i, nodeV->seq21, x->nextNode->seqi_n, routeU->vehicle);
					resultMoves[3][2] = resultMoves[2][2];
					resultMoves[3][3] = resultMoves[2][3];
				}
			}
		}
	}

	if (!nodeU->isDepot)
	{
		resultMoves[1][0] += seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq1, nodeV->seqi_n, routeV->vehicle);
		if (!x->isDepot)
		{
			resultMoves[2][0] += seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq12, nodeV->seqi_n, routeV->vehicle);
			resultMoves[3][0] += seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq21, nodeV->seqi_n, routeV->vehicle);
		}
	}

	if (!nodeV->isDepot)
	{
		if (!nodeU->isDepot)
		{
			resultMoves[1][1] += seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq1, y->seqi_n, routeV->vehicle);
			if (!x->isDepot)
			{
				resultMoves[2][1] += seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq12, y->seqi_n, routeV->vehicle);
				resultMoves[3][1] += seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq21, y->seqi_n, routeV->vehicle);
			}
		}
		if (!y->isDepot && !nodeU->isDepot && !x->isDepot)
		{
			temp = seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq12, y->nextNode->seqi_n, routeV->vehicle);
			resultMoves[2][2] += temp ;
			resultMoves[2][3] += temp ;
			temp = seq->evaluationLB(nodeV->pred->seq0_i, nodeU->seq21, y->nextNode->seqi_n, routeV->vehicle);
			resultMoves[3][2] += temp ;
			resultMoves[3][3] += temp ;
		}
	}

	/* WE IDENTIFY WHICH MOVES CAN BE FILTERED OUT USING THE LOWER BOUND */

	for (int i=0 ; i<4; i++)
	{
		for (int j=0 ; j < 4 ; j++)
		{
			shouldBeTested[i][j] = (resultMoves[i][j] < costZero) ;
		}
	}
	shouldBeTested[0][0] = true ;
	resultMoves[0][0] = costZero ;

	/* AND NOW WE TEST THE MOVES THAT HAVE A CHANCE TO BE IMPROVING */
	// Exactly the same code as previously, but using "seq->evaluation" instead of "seq->evaluationLB"

	if (!nodeU->isDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(nodeU->pred->seq0_i, x->seqi_n, routeU->vehicle);

		if (!nodeV->isDepot && shouldBeTested[1][1])
			resultMoves[1][1] = seq->evaluation(nodeU->pred->seq0_i, nodeV->seq1, x->seqi_n, routeU->vehicle);

		if (!x->isDepot)
		{
			if (shouldBeTested[2][0] || shouldBeTested[3][0]) 
			{
				resultMoves[2][0] = seq->evaluation(nodeU->pred->seq0_i, x->nextNode->seqi_n, routeU->vehicle);
				resultMoves[3][0] = resultMoves[2][0] ;
			}
			if (!nodeV->isDepot)
			{
				if (shouldBeTested[2][1] || shouldBeTested[3][1])
				{
					resultMoves[2][1] = seq->evaluation(nodeU->pred->seq0_i, nodeV->seq1, x->nextNode->seqi_n, routeU->vehicle);
					resultMoves[3][1] = resultMoves[2][1] ;
				}
				if (!y->isDepot)
				{
					if (shouldBeTested[2][2] || shouldBeTested[3][2])
					{
						resultMoves[2][2] = seq->evaluation(nodeU->pred->seq0_i, nodeV->seq12, x->nextNode->seqi_n, routeU->vehicle);
						resultMoves[3][2] = resultMoves[2][2];
					}
					if (shouldBeTested[2][3] || shouldBeTested[3][3])
					{
						resultMoves[2][3] = seq->evaluation(nodeU->pred->seq0_i, nodeV->seq21, x->nextNode->seqi_n, routeU->vehicle);
						resultMoves[3][3] = resultMoves[2][3];
					}
				}
			}
		}
	}

	if (!nodeU->isDepot)
	{
		if (shouldBeTested[1][0]) resultMoves[1][0] += seq->evaluation(nodeV->pred->seq0_i, nodeU->seq1, nodeV->seqi_n, routeV->vehicle);
		if (!x->isDepot)
		{
			if (shouldBeTested[2][0]) resultMoves[2][0] += seq->evaluation(nodeV->pred->seq0_i, nodeU->seq12, nodeV->seqi_n, routeV->vehicle);
			if (shouldBeTested[3][0]) resultMoves[3][0] += seq->evaluation(nodeV->pred->seq0_i, nodeU->seq21, nodeV->seqi_n, routeV->vehicle);
		}
	}

	if (!nodeV->isDepot)
	{
		if (!nodeU->isDepot)
		{
			if (shouldBeTested[1][1]) resultMoves[1][1] += seq->evaluation(nodeV->pred->seq0_i, nodeU->seq1, y->seqi_n, routeV->vehicle);
			if (!x->isDepot)
			{
				if (shouldBeTested[2][1]) resultMoves[2][1] += seq->evaluation(nodeV->pred->seq0_i, nodeU->seq12, y->seqi_n, routeV->vehicle);
				if (shouldBeTested[3][1]) resultMoves[3][1] += seq->evaluation(nodeV->pred->seq0_i, nodeU->seq21, y->seqi_n, routeV->vehicle);
			}
		}
		if (!y->isDepot && !nodeU->isDepot && !x->isDepot)
		{
			if (shouldBeTested[2][2] || shouldBeTested[2][3])
			{
				temp = seq->evaluation(nodeV->pred->seq0_i, nodeU->seq12, y->nextNode->seqi_n, routeV->vehicle);
				resultMoves[2][2] += temp ;
				resultMoves[2][3] += temp ;
			}
			if (shouldBeTested[3][2] || shouldBeTested[3][3])
			{
				temp = seq->evaluation(nodeV->pred->seq0_i, nodeU->seq21, y->nextNode->seqi_n, routeV->vehicle);
				resultMoves[3][2] += temp ;
				resultMoves[3][3] += temp ;
			}
		}
	}

	// We identify the best move among relocate, swap, CROSS and I-CROSS between the node pair U,V
	ibest = 0 ; jbest = 0 ;
	moveMin = 1.e30 ;
	for (int i=0 ; i<4; i++ )
	{
		for (int j=0 ; j < 4 ; j++)
		{
			if (shouldBeTested[i][j] && resultMoves[i][j] < moveMin - EPSILON_LS )
			{
				moveMin = resultMoves[i][j] ;
				ibest = i ;
				jbest = j ;
			}
		}
	}

	// If no improving move between U and V, we return
	if ( ibest == 0 && jbest == 0) 
		return 0 ;

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	Node * placeV = nodeV->pred ;
	Node * placeU = nodeU->pred ;
	reinitSingleDayMoves(placeU->route); // Say that all moves involving this route must be tested again
	reinitSingleDayMoves(placeV->route); // Say that all moves involving this route must be tested again

	// Moving the nodes in the solution structure
	if ( ibest == 1 || ibest == 2) insertNoeud(nodeU, placeV);
	if ( ibest == 2) insertNoeud(x, nodeU);
	if ( ibest == 3 ) { insertNoeud(x,placeV); insertNoeud(nodeU, x); }

	if ( jbest == 1 || jbest == 2) insertNoeud(nodeV, placeU);
	if ( jbest == 2) insertNoeud(y, nodeV);
	if ( jbest == 3 ) { insertNoeud(y,placeU); insertNoeud(nodeV, y); }

	// Update the pre-processed data on the subsequences of the route
	placeU->route->updateRouteData(false);
	placeV->route->updateRouteData(false);
    setRouteEmpty(nodeU->day); // Keep a pointer on the first empty route

	researchCompleted = false ; // Not finished the search
	nbInterSwap ++ ;
	return 1 ; // Return Success
}

int LocalSearch::interRoute2Opt ()
{
	// Testing 2-Opt* between U and V
	double cost ;
	SeqData * seq = nodeU->seq0_i ;

	if  (routeU->depot->cour != routeV->depot->cour || (nodeU->pred->isDepot && nodeV->isDepot))
		return 0 ; // Cannot do a 2-Opt* if the depot is placed in the wrong way

	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ; // Cost before the move

	// Lower bound on the move value
	cost = seq->evaluationLB(nodeU->pred->seq0_i, y->seqi_n, routeU->vehicle) + seq->evaluationLB(nodeV->seq0_i, nodeU->seqi_n, routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; // Exit if no chance of improvement

	// Exact move evaluation
	cost = seq->evaluation(nodeU->pred->seq0_i, y->seqi_n, routeU->vehicle) + seq->evaluation(nodeV->seq0_i, nodeU->seqi_n, routeV->vehicle) - costZero ;
	if ( cost  > -EPSILON_LS ) 
		return 0 ; 

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitSingleDayMoves(routeU);  // Say that all moves involving this route must be tested again
	reinitSingleDayMoves(routeV);  // Say that all moves involving this route must be tested again
	Node * tempU = nodeU ;
    nodeU = nodeU->pred ;
	x = nodeU->nextNode ;

	// Updating the solution
	Node * count ;
	Node * depotU = routeU->depot ;
	Node * depotV = routeV->depot ;
	Node * depotUFin = depotU->pred ;
	Node * depotVFin = depotV->pred ;
	Node * depotUpred = depotUFin->pred ;

	count = y ;
	while ( !count->isDepot )
	{
		count->route = routeU ;
		count = count->nextNode ;
	}

	count = x ;
	while ( !count->isDepot )
	{
		count->route = routeV ;
		count = count->nextNode ;
	}

    nodeU->nextNode = y ;
	y->pred = nodeU ;
    nodeV->nextNode = x ;
	x->pred = nodeV ;

	if (x->isDepot)
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->nextNode = depotUFin ;
        nodeV->nextNode = depotVFin ;
		depotVFin->pred = nodeV ;
	}
	else
	{
		depotUFin->pred = depotVFin->pred ;
		depotUFin->pred->nextNode = depotUFin ;
		depotVFin->pred = depotUpred ;
		depotVFin->pred->nextNode = depotVFin ;
	}

	// Update the pre-processed data on the subsequences of the route
	routeU->updateRouteData(false);
	routeV->updateRouteData(false);
    setRouteEmpty(nodeU->day); // Keep a pointer on the first empty route

	researchCompleted = false ; // Not finished the search
	nbInter2Opt ++ ;
    nodeU = tempU ;
	x = nodeU->nextNode ;
	return 1 ; // Return Success
}


int LocalSearch::interRoute2OptInv()
{
	// 2-Opt* with route inversions
	SeqData * seq = nodeU->seq0_i ;
	double cost ;
	double costTemp ;
	double costTempReverse = 1.e30 ;
	bool reverseRouteU, reverseRouteV ; 
	if (nodeU->route == nodeV->route) { return 0 ; }
	double costZero = routeU->currentRouteCost + routeV->currentRouteCost ;

	// Compute a lower bound on the value of the move
	cost = - costZero ;
	costTemp = seq->evaluationLB(y->seqn_i,x->seqi_n,routeU->vehicle) ;
	costTempReverse = seq->evaluationLB(x->seqn_i,y->seqi_n,routeU->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteU = false ; }
	else  { cost += costTempReverse ; reverseRouteU = true ; }

	costTemp = seq->evaluationLB(nodeV->seq0_i, nodeU->seqi_0, routeV->vehicle) ;
	costTempReverse = seq->evaluationLB(nodeU->seq0_i, nodeV->seqi_0, routeV->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteV = false ; }
	else  { cost += costTempReverse ; reverseRouteV = true ; }

	if ( cost  > -EPSILON_LS ) 
		return 0 ;  // Exit if no chance of improvement

	// Test the real move cost
	cost = - costZero ;
	costTemp = seq->evaluation(y->seqn_i,x->seqi_n,routeU->vehicle) ;
	costTempReverse = seq->evaluation(x->seqn_i,y->seqi_n,routeU->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteU = false ; }
	else  { cost += costTempReverse ; reverseRouteU = true ; }

	costTemp = seq->evaluation(nodeV->seq0_i, nodeU->seqi_0, routeV->vehicle) ;
	costTempReverse = seq->evaluation(nodeU->seq0_i, nodeV->seqi_0, routeV->vehicle) ;
	if (costTemp < costTempReverse) { cost += costTemp ; reverseRouteV = false ; }
	else  { cost += costTempReverse ; reverseRouteV = true ; }

	if ( cost  > -EPSILON_LS ) 
		return 0 ;

	/* AN IMPROVING MOVE HAS BEEN DETECTED, IT IS DIRECTLY APPLIED */

	reinitSingleDayMoves(nodeU->route);  // Say that all moves involving this route must be tested again
	reinitSingleDayMoves(nodeV->route);  // Say that all moves involving this route must be tested again
	Node * depotU = routeU->depot ;
	Node * depotV = routeV->depot ;
	Node * depotUFin = routeU->depot->pred ;
	Node * depotVFin = routeV->depot->pred ;
	Node * depotUSuiv = depotU->nextNode ;

	// Update the solution
	Node * temp ;
	Node * yy = y ;
	Node * uu = nodeU ;

	while ( !yy->isDepot )
	{
		temp = yy->nextNode ;
		yy->nextNode = yy->pred ;
		yy->pred = temp ;
		yy->route = routeU ;
		yy = temp ;
	}

	while ( !uu->isDepot )
	{
		temp = uu->pred ;
		uu->pred = uu->nextNode ;
		uu->nextNode = temp ;
		uu->route = routeV ;
		uu = temp ;
	}

    nodeV->nextNode = nodeU ;
    nodeU->pred = nodeV ;
	y->nextNode = x ;
	x->pred = y ;

	if (y->isDepot)
	{
		depotVFin->nextNode = depotV ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->nextNode = depotVFin ;
		depotU->nextNode = x ;
		x->pred = depotU ;
	}
	else if (nodeU->isDepot )
	{
		depotU->nextNode = depotVFin->pred ;
		depotU->nextNode->pred = depotU ;
		depotU->pred = depotUFin ;
		depotVFin->pred = nodeV ;
        nodeV->nextNode = depotVFin ;
	}
	else
	{
		depotU->nextNode = depotVFin->pred ;
		depotU->nextNode->pred = depotU ;
		depotVFin->pred = depotUSuiv ;
		depotVFin->pred->nextNode = depotVFin ;
	}

	// Reverse if needed
	if (reverseRouteU) routeU->reverse();
	if (reverseRouteV) routeV->reverse();

	// Update the pre-processed data on the subsequences of the route
	routeU->updateRouteData(false);
	routeV->updateRouteData(false);
    setRouteEmpty(nodeU->day); // Keep a pointer on the first empty route

	researchCompleted = false ; // Not finished the search
	nbInter2Opt ++ ;
	return 1 ; // Return Success
}

int LocalSearch::intraRouteGeneralInsertDroite ()
{
	// For a pair of nodes U, V, IN THE SAME ROUTE, tests together the Relocate, Swap, and variants of CROSS and I-CROSS limited to two consecutive nodes.
	Node * tempU = nodeU ;
	Node * tempV = nodeV ;
	bool turned = false ;

	int decalage = nodeV->place - nodeU->place ; // Decalage is the difference of index between U and V
	if (decalage >= -1 && decalage <= 1) return 0 ; // This means that U and V are already consecutive, testing these moves is useless

	// If decalage < 0, this means that V is on the left of U in the route
	// We don't want to deal with this case, so we simply inverse the notations of U and V 
	// And we deal in the following with only the case where V is on the right
	if ( decalage < 0 )
	{
        nodeU = tempV ;
        nodeV = tempU ;
		x = nodeU->nextNode ;
		y = nodeV->nextNode ;
		decalage = -decalage ;
		turned = true ;
	}

	double moveMin ;
	int ibest, jbest ;
	SeqData * seq = nodeU->seq0_i ;
	Node * Upred = nodeU->pred ;
	Node * placeU = nodeU->pred ;
	Node * placeV = nodeV->pred ;
	Node * xsuiv = x->nextNode ;
	Node * ysuiv = y->nextNode ;
	double costZero = routeU->currentRouteCost ;

	// This table will keep the result of the moves
	// 0 -> send nothing
	// 1 -> send U
	// 2 -> send U,Unext
	// 3 -> send Unext,U
	for (int i=0 ; i<4 ; i++ )
		for (int j=0 ; j<4 ; j++ )
			resultMoves[i][j] = 1.e29 ;

	resultMoves[0][0] = costZero ;

	if (decalage >= 3) // General case
	{
		if (!nodeV->isDepot && turned) // relocate of V before U (V cannot be a depot, but U can be)
		{
			mySeqs.clear();
			mySeqs.push_back(Upred->seq0_i);
			mySeqs.push_back(nodeV->seq1);
			// The intra-route moves are the only moves which can involve general subsequences (i,j) in the middle of the route
			// To reduce a bit the pre-processing, we keep only sequences of limited size (<= 15), even if this involves to concatenate more pieces
			// The concatenation of the good number of pieces is done by "addSeqDataInPieces"
			addSeqDataInPieces(nodeU, nodeV->place - 1 - nodeU->place, nodeU->day);
			mySeqs.push_back(y->seqi_n);
			resultMoves[0][1] = seq->evaluationLB(mySeqs, routeU->vehicle);
			if (!y->isDepot) // exchange U and (V,Y or Y,V)
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				mySeqs.push_back(nodeV->seq12);
				addSeqDataInPieces(nodeU, nodeV->place - 1 - nodeU->place, nodeU->day);
				mySeqs.push_back(ysuiv->seqi_n);
				resultMoves[0][2] = seq->evaluationLB(mySeqs, routeU->vehicle);
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				mySeqs.push_back(nodeV->seq21);
				addSeqDataInPieces(nodeU, nodeV->place - 1 - nodeU->place, nodeU->day);
				mySeqs.push_back(ysuiv->seqi_n);
				resultMoves[0][3] = seq->evaluationLB(mySeqs, routeU->vehicle);
			}
		}

		if (!turned)
		{
			mySeqs.clear();
			mySeqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
			mySeqs.push_back(nodeU->seq1);
			mySeqs.push_back(nodeV->seqi_n);
			resultMoves[1][0] = seq->evaluationLB(mySeqs, routeU->vehicle);
		}
		if (!nodeV->isDepot) // exchange U and V
		{
			mySeqs.clear();
			mySeqs.push_back(Upred->seq0_i);
			mySeqs.push_back(nodeV->seq1);
			addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
			mySeqs.push_back(nodeU->seq1);
			mySeqs.push_back(y->seqi_n);
			resultMoves[1][1] = seq->evaluationLB(mySeqs, routeU->vehicle);
			if (!y->isDepot && turned) // exchange U and (V,Y or Y,V)
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				mySeqs.push_back(nodeV->seq12);
				addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
				mySeqs.push_back(nodeU->seq1);
				mySeqs.push_back(ysuiv->seqi_n);
				resultMoves[1][2] = seq->evaluationLB(mySeqs, routeU->vehicle);

				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				mySeqs.push_back(nodeV->seq21);
				addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
				mySeqs.push_back(nodeU->seq1);
				mySeqs.push_back(ysuiv->seqi_n);
				resultMoves[1][3] = seq->evaluationLB(mySeqs, routeU->vehicle);
			}
		}

		if (!x->isDepot)
		{
			if (!turned)
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
				mySeqs.push_back(nodeU->seq12);
				mySeqs.push_back(nodeV->seqi_n);
				resultMoves[2][0] = seq->evaluationLB(mySeqs, routeU->vehicle);
			}

			if (!nodeV->isDepot) // exchange U,X and V
			{
				if (!turned)
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq1);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq12);
					mySeqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluationLB(mySeqs, routeU->vehicle);
				}
				if (!y->isDepot) // exchange U,X and (V,Y or Y,V)
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq12);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq12);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[2][2] = seq->evaluationLB(mySeqs, routeU->vehicle);

					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq21);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq12);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[2][3] = seq->evaluationLB(mySeqs, routeU->vehicle);
				}
			}

			if (!turned)
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
				mySeqs.push_back(nodeU->seq21);
				mySeqs.push_back(nodeV->seqi_n);
				resultMoves[3][0] = seq->evaluationLB(mySeqs, routeU->vehicle);
			}
			if (!nodeV->isDepot) // exchange X,U and V
			{
				if (!turned)
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq1);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq21);
					mySeqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluationLB(mySeqs, routeU->vehicle);
				}
				if (!y->isDepot) // exchange X,U and (V,Y or Y,V)
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq12);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq21);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[3][2] = seq->evaluationLB(mySeqs, routeU->vehicle);

					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq21);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq21);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[3][3] = seq->evaluationLB(mySeqs, routeU->vehicle);
				}
			}
		}
	}
	else if (decalage == 2) // U and V are almost consecutive, we are in this situation : UXVY. Useful moves must be tested as special cases
	{
		// XUVY
		resultMoves[1][0] = seq->evaluationLB(Upred->seq0_i, nodeU->seq21, nodeV->seqi_n, routeU->vehicle);
		if (!nodeV->isDepot)
		{
			// VXUY
			resultMoves[1][1] = seq->evaluationLB(Upred->seq0_i, nodeV->seq1, nodeU->seq21, y->seqi_n, routeU->vehicle);
			// VUXY
			resultMoves[2][1] = seq->evaluationLB(Upred->seq0_i, nodeV->seq1, nodeU->seq12, y->seqi_n, routeU->vehicle);
			if (!y->isDepot)
			{
				// VYXU
				resultMoves[1][2] = seq->evaluationLB(Upred->seq0_i, nodeV->seq12, nodeU->seq21, ysuiv->seqi_n, routeU->vehicle);
				// YVXU
				resultMoves[1][3] = seq->evaluationLB(Upred->seq0_i, nodeV->seq21, nodeU->seq21, ysuiv->seqi_n, routeU->vehicle);
				// VYUX
				resultMoves[2][2] = seq->evaluationLB(Upred->seq0_i, nodeV->seq12, nodeU->seq12, ysuiv->seqi_n, routeU->vehicle);
				// YVUX
				resultMoves[2][3] = seq->evaluationLB(Upred->seq0_i, nodeV->seq21, nodeU->seq12, ysuiv->seqi_n, routeU->vehicle);
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// HERE STARTS THE SECOND PHASE AFTER PREPROCESSING
	// ONLY NON-FILTERED MOVES THAT HAVE A CHANCE TO BE IMPROVING
	/////////////////////////////////////////////////////////////

	for (int i=0 ; i<4 ; i++ )
	{
		for (int j=0 ; j<4 ; j++ )
		{
			shouldBeTested[i][j] = (resultMoves[i][j] < costZero) ;
		}
	}
	shouldBeTested[0][0] = true  ;
	resultMoves[0][0] = costZero ;


	// Same procedure as previously, but with "seq->evaluation" instead of "seq->evaluationLB"
	if (decalage >= 3)
	{
		if (!nodeV->isDepot && turned)
		{

			if (shouldBeTested[0][1])
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				mySeqs.push_back(nodeV->seq1);
				addSeqDataInPieces(nodeU, nodeV->place - 1 - nodeU->place, nodeU->day);
				mySeqs.push_back(y->seqi_n);
				resultMoves[0][1] = seq->evaluation(mySeqs, routeU->vehicle);
			}

			if (!y->isDepot) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[0][2])
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq12);
					addSeqDataInPieces(nodeU, nodeV->place - 1 - nodeU->place, nodeU->day);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[0][2] = seq->evaluation(mySeqs, routeU->vehicle);
				}

				if (shouldBeTested[0][3])
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq21);
					addSeqDataInPieces(nodeU, nodeV->place - 1 - nodeU->place, nodeU->day);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[0][3] = seq->evaluation(mySeqs, routeU->vehicle);
				}
			}
		}

		if (!turned && shouldBeTested[1][0])
		{
			mySeqs.clear();
			mySeqs.push_back(Upred->seq0_i);
			addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
			mySeqs.push_back(nodeU->seq1);
			mySeqs.push_back(nodeV->seqi_n);
			resultMoves[1][0] = seq->evaluation(mySeqs, routeU->vehicle);
		}
		if (!nodeV->isDepot) // exchange U and V
		{
			if (shouldBeTested[1][1])
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				mySeqs.push_back(nodeV->seq1);
				addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
				mySeqs.push_back(nodeU->seq1);
				mySeqs.push_back(y->seqi_n);
				resultMoves[1][1] = seq->evaluation(mySeqs, routeU->vehicle);
			}

			if (!y->isDepot && turned) // exchange U and (V,Y or Y,V)
			{
				if (shouldBeTested[1][2])
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq12);
					addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
					mySeqs.push_back(nodeU->seq1);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[1][2] = seq->evaluation(mySeqs, routeU->vehicle);
				}

				if (shouldBeTested[1][3])
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq21);
					addSeqDataInPieces(x, nodeV->place - 1 - x->place, nodeU->day);
					mySeqs.push_back(nodeU->seq1);
					mySeqs.push_back(ysuiv->seqi_n);
					resultMoves[1][3] = seq->evaluation(mySeqs, routeU->vehicle);
				}
			}
		}

		if (!x->isDepot)
		{
			if (shouldBeTested[2][0] && !turned)
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
				mySeqs.push_back(nodeU->seq12);
				mySeqs.push_back(nodeV->seqi_n);
				resultMoves[2][0] = seq->evaluation(mySeqs, routeU->vehicle);
			}

			if (!nodeV->isDepot) // exchange U,X and V
			{
				if (shouldBeTested[2][1] && !turned)
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq1);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq12);
					mySeqs.push_back(y->seqi_n);
					resultMoves[2][1] = seq->evaluation(mySeqs, routeU->vehicle);
				}
				if (!y->isDepot) // exchange U,X and (V,Y or Y,V)
				{
					if (shouldBeTested[2][2])
					{
						mySeqs.clear();
						mySeqs.push_back(Upred->seq0_i);
						mySeqs.push_back(nodeV->seq12);
						addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
						mySeqs.push_back(nodeU->seq12);
						mySeqs.push_back(ysuiv->seqi_n);
						resultMoves[2][2] = seq->evaluation(mySeqs, routeU->vehicle);
					}

					if (shouldBeTested[2][3])
					{
						mySeqs.clear();
						mySeqs.push_back(Upred->seq0_i);
						mySeqs.push_back(nodeV->seq21);
						addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
						mySeqs.push_back(nodeU->seq12);
						mySeqs.push_back(ysuiv->seqi_n);
						resultMoves[2][3] = seq->evaluation(mySeqs, routeU->vehicle);
					}
				}
			}

			if (shouldBeTested[3][0] && !turned)
			{
				mySeqs.clear();
				mySeqs.push_back(Upred->seq0_i);
				addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
				mySeqs.push_back(nodeU->seq21);
				mySeqs.push_back(nodeV->seqi_n);
				resultMoves[3][0] = seq->evaluation(mySeqs, routeU->vehicle);
			}
			if (!nodeV->isDepot) // exchange X,U and V
			{
				if (shouldBeTested[3][1] && !turned)
				{
					mySeqs.clear();
					mySeqs.push_back(Upred->seq0_i);
					mySeqs.push_back(nodeV->seq1);
					addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
					mySeqs.push_back(nodeU->seq21);
					mySeqs.push_back(y->seqi_n);
					resultMoves[3][1] = seq->evaluation(mySeqs, routeU->vehicle);
				}
				if (!y->isDepot) // exchange X,U and (V,Y or Y,V)
				{
					if (shouldBeTested[3][2])
					{
						mySeqs.clear();
						mySeqs.push_back(Upred->seq0_i);
						mySeqs.push_back(nodeV->seq12);
						addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
						mySeqs.push_back(nodeU->seq21);
						mySeqs.push_back(ysuiv->seqi_n);
						resultMoves[3][2] = seq->evaluation(mySeqs, routeU->vehicle);
					}

					if (shouldBeTested[3][3])
					{
						mySeqs.clear();
						mySeqs.push_back(Upred->seq0_i);
						mySeqs.push_back(nodeV->seq21);
						addSeqDataInPieces(xsuiv, nodeV->place - 1 - xsuiv->place, nodeU->day);
						mySeqs.push_back(nodeU->seq21);
						mySeqs.push_back(ysuiv->seqi_n);
						resultMoves[3][3] = seq->evaluation(mySeqs, routeU->vehicle);
					}
				}
			}
		}
	}
	else if (decalage == 2)
	{
		// XUVY
		if (shouldBeTested[1][0]) resultMoves[1][0] = seq->evaluation(Upred->seq0_i, nodeU->seq21, nodeV->seqi_n, routeU->vehicle);
		if (!nodeV->isDepot)
		{
			// VXUY
			if (shouldBeTested[1][1]) resultMoves[1][1] = seq->evaluation(Upred->seq0_i, nodeV->seq1, nodeU->seq21, y->seqi_n, routeU->vehicle);
			// VUXY
			if (shouldBeTested[2][1]) resultMoves[2][1] = seq->evaluation(Upred->seq0_i, nodeV->seq1, nodeU->seq12, y->seqi_n, routeU->vehicle);
			if (!y->isDepot)
			{
				// VYXU
				if (shouldBeTested[1][2]) resultMoves[1][2] = seq->evaluation(Upred->seq0_i, nodeV->seq12, nodeU->seq21, ysuiv->seqi_n, routeU->vehicle);
				// YVXU
				if (shouldBeTested[1][3]) resultMoves[1][3] = seq->evaluation(Upred->seq0_i, nodeV->seq21, nodeU->seq21, ysuiv->seqi_n, routeU->vehicle);
				// VYUX
				if (shouldBeTested[2][2]) resultMoves[2][2] = seq->evaluation(Upred->seq0_i, nodeV->seq12, nodeU->seq12, ysuiv->seqi_n, routeU->vehicle);
				// YVUX
				if (shouldBeTested[2][3]) resultMoves[2][3] = seq->evaluation(Upred->seq0_i, nodeV->seq21, nodeU->seq12, ysuiv->seqi_n, routeU->vehicle);
			}
		}
	}

	// Identify the best move involving U,V and apply
	ibest = 0 ;
	jbest = 0 ;
	moveMin = 1.e30 ;
	for (int i=0 ; i<4 ; i++)
	{
		for (int j=0 ; j<4 ; j++ )
		{
			if (shouldBeTested[i][j] && resultMoves[i][j] < moveMin - EPSILON_LS )
			{
				moveMin = resultMoves[i][j] ;
				ibest = i ;
				jbest = j ;
			}
		}
	}

	// No improving move has been found
	if ( ibest == 0 && jbest == 0 ) 
	{
        nodeU = tempU ;
        nodeV = tempV ;
		x = nodeU->nextNode ;
		y = nodeV->nextNode ;
		return 0 ;
	}


	///////////////////////////////////////////
	// HERE STARTS THE MOVE APPLICATION PROCESS
	///////////////////////////////////////////

	reinitSingleDayMoves(routeU);  // Say that all moves involving this route must be tested again

	// Update the solution
	if (decalage >= 3)
	{
		if ( ibest == 1 || ibest == 2) insertNoeud(nodeU, placeV);
		if ( ibest == 2 ) insertNoeud(x, nodeU);
		if ( ibest == 3 ) { insertNoeud(x,placeV); insertNoeud(nodeU, x); }

		if ( jbest == 1 || jbest == 2) insertNoeud(nodeV, placeU);
		if ( jbest == 2) insertNoeud(y, nodeV);
		if ( jbest == 3 ) { insertNoeud(y,placeU); insertNoeud(nodeV, y); }
	}

	// Special cases of decalage == 2
	else if (decalage == 2)
	{
		if (ibest == 1 && jbest == 0) { insertNoeud(nodeU, x); }
		else if (ibest == 1 && jbest == 1) { insertNoeud(x, nodeV); insertNoeud(nodeU, x);  }
		else if (ibest == 1 && jbest == 2) { insertNoeud(x,y); insertNoeud(nodeU, x);  }
		else if (ibest == 1 && jbest == 3) { insertNoeud(nodeV, y); insertNoeud(x, nodeV); insertNoeud(nodeU, x);  }
		else if (ibest == 2 && jbest == 1) { insertNoeud(nodeV, nodeU->pred) ;  }
		else if (ibest == 2 && jbest == 2) { insertNoeud(nodeV, nodeU->pred) ;  insertNoeud(y, nodeV) ; }
		else if (ibest == 2 && jbest == 3) { insertNoeud(y, nodeU->pred) ;  insertNoeud(nodeV, y) ; }
		else throw string ("ERROR move intra-route") ;
	} 

	routeU->updateRouteData(false); // Update the pre-processed data on the subsequences of the route
    setRouteEmpty(nodeU->day); // Keep a pointer on the first empty route
	researchCompleted = false ; // Not finished the search
	nbIntraSwap ++ ;
    nodeU = tempU ;
    nodeV = tempV ;
	x = nodeU->nextNode ;
	y = nodeV->nextNode ;
	return 1 ; // Return Success
}


int LocalSearch::intraRoute2Opt ()
{
	// Evaluation procedure for 2-Opt
	Node * nodeNum = nodeU->nextNode ;
	Node * nodeUpred = nodeU->pred ;
	Node * temp ;
	SeqData * seq = nodeU->seq0_i ;

	double cost ;
	double costZero = routeU->currentRouteCost ;

	mySeqs.clear();
	mySeqs.push_back(nodeU->pred->seq0_i);
	addReverseSeqDataInPieces(nodeU, nodeV->place - nodeU->place, nodeU->day);
	mySeqs.push_back(nodeV->nextNode->seqi_n);

	// Compute the lower bound on move value and exits if no possible improvement
	cost = seq->evaluationLB(mySeqs, routeU->vehicle) ;
	if (cost - costZero  > -EPSILON_LS)  
		return 0 ;

	// Compute the real move value
	cost = seq->evaluation(mySeqs, routeU->vehicle) ;
	if (cost - costZero  > -EPSILON_LS)  
		return 0 ;

	// Apply the move and update the solution
	reinitSingleDayMoves(routeU);
    nodeU->pred = nodeNum ;
    nodeU->nextNode = y ;

	while (nodeNum != nodeV )
	{
		temp = nodeNum->nextNode ;
		nodeNum->nextNode = nodeNum->pred ;
		nodeNum->pred = temp ;
		nodeNum = temp ;
	}

    nodeV->nextNode = nodeV->pred ;
    nodeV->pred = nodeUpred ;
	nodeUpred->nextNode = nodeV ;
	y->pred = nodeU ;
	routeU->updateRouteData(false);  // Update the pre-processed data on the subsequences of the route
    setRouteEmpty(nodeU->day); // Keep a pointer on the first empty route
	researchCompleted = false ; // Not finished the search
	nbIntra2Opt ++ ;
	return 1 ; // Return Success
}

int LocalSearch::searchBetterPattern (int client)
{
	pattern pattern1 = individu->chromP[client] ;
	pattern pattern2, meilleurPattern ;
	int indexMeilleur = -1 ;
	meilleurPattern.pat = -1000000 ;
	int temp, calcul, depot ;
	double depense = 0 ;
	double meilleureDepense = 1.e30 ;
	Node * noeudTravail ;
    intraDayDisplacement = false ; 	// this flag is only raised in case there is a better insertion in the same day

	for (int pat = 0 ; pat < (int)params->cli[client].visits.size() ; pat ++)
	{
		// testing a new pattern
		pattern2 = params->cli[client].visits[pat] ;
		testingIncumbentPattern = (individu->chromP[client].pat == pattern2.pat && individu->chromP[client].dep == pattern2.dep) ;
		calcul = pattern2.pat ;
		depot = pattern2.dep ;
		depense = pattern2.cost ;
		for (int k = 0 ; k < params->formerNbDays ; k++)
		{
			temp = calcul % 2 ;
			calcul = calcul/2 ;
			if (temp == 1)
			{
				noeudTravail = &clients[params->formerNbDays - k + depot * params->formerNbDays][client] ;
				// If the insertion in this day has not been computed yet 
				// (Note that it would be possible to do way faster by considering the fact that only the demand of 
				// a single delivery changes... and sometimes its even the same for different patterns)
				// Still, this is a simplistic implementation for the tests PCARP, no need for the highest performance
				if (noeudTravail->costInsertion[pat] > 1.e29 || firstLoop)
                    computeCostInsertion(noeudTravail, pat) ;
				depense += noeudTravail->costInsertion[pat] ;
			}
		}

		if (depense < meilleureDepense - EPSILON_LS
			|| (depense < meilleureDepense + EPSILON_LS && pattern2.dep == pattern1.dep && pattern2.pat == pattern1.pat) )
		{
			meilleureDepense = depense ;
			meilleurPattern = pattern2 ;
			indexMeilleur = pat ;
		}
	}

	if (meilleurPattern.pat == -1000000) 
		throw string ("ERROR when computing the best pattern !") ;

	// Applying the move if a better pattern has been found
	if ( meilleurPattern.pat != pattern1.pat || meilleurPattern.dep != pattern1.dep || intraDayDisplacement)
	{
		// removing the current occurences of this customer
		calcul = pattern1.pat ;
		for (int k = 0 ; k < params->formerNbDays ; k++)
		{
			temp = calcul % 2 ;
			calcul = calcul/2 ;
			if (temp == 1)
				removeNoeud(&clients[params->formerNbDays - k + pattern1.dep * params->formerNbDays][client]);
		}

		// (PCARP) Updating the chromP (necessary to do now, not later, otherwise the wrong data is set to pre-process the SeqData)
		// When adding the nodes in the next block of instructions
		individu->chromP[client] = meilleurPattern ;

		// Inserting in the new locations
		calcul = meilleurPattern.pat ;
		depot = meilleurPattern.dep ;
		for (int k = 0 ; k < params->formerNbDays ; k++)
		{
			temp = calcul % 2 ;
			calcul = calcul/2 ;
			if (temp == 1)
			{
				Node * hereCli = &clients[params->formerNbDays - k + meilleurPattern.dep * params->formerNbDays][client] ;
				Node * thereCli = hereCli->placeInsertion[indexMeilleur] ;
				addNoeud(hereCli,thereCli);
			}	
		}
		//cout << "Inserting Node " << client << " With pattern " << meilleurPattern.pat << " Flag is : " << intraDayDisplacement << endl ;

		researchCompleted = false ;
		return 1 ;
	}
	else return 0 ;
}

void LocalSearch::computeCostInsertion(Node * client, int pattern)
{
	// Computing the best cost for inserting a client in its day
	Route * myRoute ;
	client->costInsertion[pattern] = 1.e30 ;
	client->placeInsertion[pattern] = NULL ;

    nodeU = client ;
	x = nodeU->nextNode ;
	routeU = nodeU->route ;

	// Find the best insertion for each route
	for (int r=0 ; r < params->numberVehicle[client->day] ; r++)
	{
		myRoute = &routes[client->day][r] ;
		if ( myRoute->coutInsertionClient[client->cour][pattern] > 1.e29 || firstLoop ) 
			evalInsertClient(myRoute,client,pattern) ;

		if ( myRoute->coutInsertionClient[client->cour][pattern] < client->costInsertion[pattern] - EPSILON_LS)
		{
			client->costInsertion[pattern] = myRoute->coutInsertionClient[client->cour][pattern] ;
			client->placeInsertion[pattern] = myRoute->placeInsertionClient[client->cour][pattern] ;
		}
	}

	// If its possible to improve the placement of a customer in a day where its already placed, and according to its current pattern, then we raise this flag
	if (client->isPresent // its placed here
		&& testingIncumbentPattern // with the same pattern
		&& client->route->coutInsertionClient[client->cour][pattern] > client->costInsertion[pattern] + EPSILON_LS) // but we can do better
		intraDayDisplacement = true ;
}

void LocalSearch::evalInsertClient (Route * R, Node * U, int pattern)
{
	// Computing the least cost for inserting a client U in a route R
	SeqData * seq = U->seq0_i ;
	Node * courNoeud ;
	double leastCost = 1.e30 ;
	double cost ;
	bool firstLoopDep = true ;

	// Tweaking the code to work with the PCARP
	// We need to account for the fact that the demand may change based on the pattern choice
	// Here doing something a bit ugly, which is to store and modify from outside the pre-processed Seqdata "U->seq1"
	// to include the good delivery quantity
	double tempDemand = U->seq1->load ;
	U->seq1->load = params->cli[U->cour].demandPatDay[params->cli[U->cour].visits[pattern].pat][U->day];

	// Some memory structures to avoid recomputing these things again and again
	R->coutInsertionClient[U->cour][pattern] = 1.e30 ;
	R->placeInsertionClient[U->cour][pattern] = NULL ;

	if (!(U->route == R) || !U->isPresent)
	{
		// Case 1 : U is not in the route
		courNoeud = R->depot->nextNode ;
		while (!courNoeud->isDepot || firstLoopDep)
		{
			if (courNoeud->isDepot) firstLoopDep = false ;
			cost = seq->evaluation(courNoeud->pred->seq0_i,U->seq1,courNoeud->seqi_n,R->vehicle);
			if ( cost < leastCost )
			{
				leastCost = cost ;
				R->placeInsertionClient[U->cour][pattern] = courNoeud->pred ;
			}
			courNoeud = courNoeud->nextNode ;
		}
		R->coutInsertionClient[U->cour][pattern] = leastCost - seq->evaluation(R->depot->pred->seq0_i,R->vehicle);
	}
	else
	{
		// Case 2 : U is already in the route R
		leastCost = seq->evaluation(U->pred->seq0_i, U->seq1, U->nextNode->seqi_n, R->vehicle);
		R->placeInsertionClient[U->cour][pattern] = U->pred ;
		courNoeud = R->depot->nextNode ;
		while (!courNoeud->isDepot || firstLoopDep)
		{
			if (courNoeud->isDepot) firstLoopDep = false ;

			if (courNoeud->place < U->place)
			{
				mySeqs.clear();
				mySeqs.push_back(courNoeud->pred->seq0_i);
				mySeqs.push_back(U->seq1);
				addSeqDataInPieces(courNoeud,U->place-1-courNoeud->place,courNoeud->day);
				mySeqs.push_back(U->nextNode->seqi_n);
				cost = seq->evaluation(mySeqs, R->vehicle);
			}
			else if (courNoeud->place > U->place + 1)
			{
				mySeqs.clear();
				mySeqs.push_back(U->pred->seq0_i);
				addSeqDataInPieces(U->nextNode, courNoeud->place - 1 - U->nextNode->place, courNoeud->day);
				mySeqs.push_back(U->seq1);
				mySeqs.push_back(courNoeud->seqi_n);
				cost = seq->evaluation(mySeqs, R->vehicle);
			}
			else cost = 1.e30 ;

			if ( cost < leastCost - EPSILON_LS )
			{
				leastCost = cost ;
				R->placeInsertionClient[U->cour][pattern] = courNoeud->pred ;
			}
			courNoeud = courNoeud->nextNode ;
		}
		R->coutInsertionClient[U->cour][pattern] = leastCost - seq->evaluation(U->pred->seq0_i, U->nextNode->seqi_n, R->vehicle);
	}

	// PCARP, setting back the pre-processed demand to its correct value
	U->seq1->load = tempDemand ;
}

void LocalSearch::shuffleRouteOrder() {
    // Shuffling the routeOrder vector
    int j, temp;
    for (int k = 0; k <= params->nbDays; k++) {
        for (int i = 0; i < (int) routeOrder[k].size() - 1; i++) {
            j = i + rand() % ((int) routeOrder[k].size() - i);
            temp = routeOrder[k][i];
            routeOrder[k][i] = routeOrder[k][j];
            routeOrder[k][j] = temp;
        }
    }
}

// 总体更新的是clients[day][clientIndex]中的moves结构，在判断closestNeighbor时，将该client视为close的，并且在当天需要服务的客户节点，作为待选择的动作
void LocalSearch::updateMoves ()
{
	int client, client2, size ;
	for (int k=1 ; k<=params->nbDays ; k++) {
        for (int i = 0; i < (int) routeOrder[k].size(); i++) {
            client = routeOrder[k][i];
            clients[k][client].moves.clear();
            size = (int) params->cli[client].neighborsCloseBefore.size();
            for (int a1 = 0; a1 < size; a1++) {
                client2 = params->cli[client].neighborsCloseBefore[a1];
                if (client2 >= params->nbDepots && clients[k][client2].isPresent)
                    clients[k][client].moves.push_back(client2);
            }
        }
	}
}

void LocalSearch::setRouteEmpty(int day)
{
    /*
     * 因为在updateLS阶段，路径是从后往前分配的，如果split切出的路径数量小雨车辆数，前面的若干的路径是空路径，
     * 空路径的特点，起点depot的下一个节点也是depot，并且起点depot的route为NULL
     * 该方法使得routeEmpty[day]指向最后一个空路径
     */
    routeEmpty[day] = NULL ;
	int route = 0 ;
	while (routeEmpty[day] == NULL && route < params->nbVehiclesPerDepot )
	{
		if (depots[day][route].nextNode->isDepot) routeEmpty[day] = depots[day][route].route ;
		route ++ ;
	}
}

void LocalSearch::nodeTestedForEachRoute (int cli, int day)
{
	for (int route = 0 ; route < (int)params->numberVehicle[day] ; route ++)
		routes[day][route].nodeAndRouteTested[cli]=true ;	
}

void LocalSearch::placeMissing ()
{
	// Insertion procedure for the crossover PIX
	int k ;
	pattern betterPattern ;
	double depense, betterDepense ;
	int indexBetter ;
	Node * nodeTravail ;
	int calcul1, calcul2 ;
	pattern pattern1, pattern2 ;

	firstLoop = true ;

	for (int day = 1 ; day <= params->nbDays ; day++)
		for (int r=0 ; r < params->numberVehicle[day] ; r++)
			routes[day][r].initiateInsertions() ;

	// We iterate on missing visits
	for (int i=0 ; i < (int)individu->toPlace.size() ; i++ )
	{
        betterDepense = 1.e30 ;
		k = individu->toPlace[i] ;
		pattern1 = individu->chromP[k] ;

		for (int pat = 0 ; pat < (int)params->cli[k].visits.size() ; pat ++)
		{
			// testing a new pattern
			pattern2 = params->cli[k].visits[pat] ;
			if (pattern::isSubset(pattern1,pattern2))
			{
				calcul2 = pattern2.pat ;
				depense = pattern2.cost ;
				for (int kk = 0 ; kk < params->formerNbDays ; kk++)
				{
					if (calcul2 % 2 == 1)
					{
                        nodeTravail = &clients[params->formerNbDays - kk + pattern2.dep * params->formerNbDays][k] ;
                        computeCostInsertion(nodeTravail, pat) ;
						depense += nodeTravail->costInsertion[pat] ;
					}
					calcul2 = calcul2/2 ;
				}

				if (depense < betterDepense)
				{
                    betterDepense = depense ;
                    betterPattern = pattern2 ;
                    indexBetter = pat ;
				}
			}
		}

		// Updating the chromP with the chosen pattern
		individu->chromP[k] = betterPattern ;

		// Inserting in the new locations
		calcul1 = pattern1.pat ;
		calcul2 = betterPattern.pat ;
		for (int kk = 0 ; kk < params->formerNbDays ; kk++)
		{
			if (calcul2 % 2 == 1 && calcul1 % 2 == 0)
			{
				Node * hereCli = &clients[params->formerNbDays - kk + betterPattern.dep * params->formerNbDays][k] ;
				Node * thereCli = hereCli->placeInsertion[indexBetter] ;
				addNoeud(hereCli,thereCli);
			}
			calcul1 = calcul1/2 ;
			calcul2 = calcul2/2 ;
		}
	}
}

void LocalSearch::insertNoeud(Node * U, Node * V)
{
	if (U->pred != V && U != V)
	{
		U->pred->nextNode = U->nextNode ;
		U->nextNode->pred = U->pred ;
		V->nextNode->pred = U ;
		U->pred = V ;
		U->nextNode = V->nextNode ;
		V->nextNode = U ;
		U->route = V->route ;
	}
}

void LocalSearch::removeNoeud(Node * U)
{
	Node * temp =  U->nextNode ;
	reinitSingleDayMoves(U->route);
	U->pred->nextNode = U->nextNode ;
	U->nextNode->pred = U->pred ;
	U->route = NULL ;
	U->pred = NULL ;
	U->nextNode = NULL ;
	temp->route->updateRouteData(false);
    setRouteEmpty(U->day);

	// Managing the other data structures
	individu->chromT[U->day].pop_back();
	removeOP(U->day, U->cour);
	U->isPresent = false ;

	// Say that the insertions on this day need to be computed again
	temp->route->initiateInsertions();
	for (int i= params->nbDepots ; i< params->nbDepots + params->nbClients ; i++)
		for (int p=0 ; p < (int)params->cli[i].visits.size() ; p++)
			clients[U->day][i].costInsertion[p] = 1.e30 ;
}

void LocalSearch::addNoeud(Node * U, Node * V)
{
	reinitSingleDayMoves(V->route);
	V->nextNode->pred = U ;
	U->pred = V;
	U->nextNode = V->nextNode ;
	V->nextNode = U ;

	// Update the routes
	U->isPresent = true ;
	U->route = V->route ;
	U->route->updateRouteData(false);
    setRouteEmpty(U->day);

	// Manage the other data structures
	individu->chromT[U->day].push_back(0);
	addOP(U->day, U->cour);

	// Say that the insertions on this day need to be computed again
	U->route->initiateInsertions();
	for (int i= params->nbDepots ; i< params->nbDepots + params->nbClients ; i++)
		for (int p=0 ; p < (int)params->cli[i].visits.size() ; p++)
			clients[U->day][i].costInsertion[p] = 1.e30 ;
}

void LocalSearch::removeOP (int day, int client) {
    // Remove one occurence in the routeOrder structure
    int it = 0;
    while (routeOrder[day][it] != client) { it++; }
    routeOrder[day][it] = routeOrder[day][(int) routeOrder[day].size() - 1];
    routeOrder[day].pop_back();
}

void LocalSearch::addOP (int day, int client) {
    // Add one element in the routeOrder structure
    int it, temp2;
    if (routeOrder[day].size() != 0) {
        it = (int) rand() % routeOrder[day].size();
        temp2 = routeOrder[day][it];
        routeOrder[day][it] = client;
        routeOrder[day].push_back(temp2);
    } else
        routeOrder[day].push_back(client);
}

void LocalSearch::controlRoutes ()
{
	bool tracesDebug = false; // set to true to activate traces
	Node * myNoeud ;
	Node * autreNoeud ;
	bool firstIt = true ;
	for (int k=1 ; k<=params->nbDays ; k++)
		for (int r=0 ; r < params->numberVehicle[k] ; r++)
		{
			firstIt = true ;
			if (tracesDebug) cout << "day(" << k << ") " << "route(" << r << ") // " ;
			myNoeud = routes[k][r].depot ;
			if (tracesDebug) cout << myNoeud->cour ;
			while (!myNoeud->isDepot || firstIt )
			{
				firstIt = false ;
				autreNoeud = myNoeud ;
				myNoeud = myNoeud->nextNode ;
				if (tracesDebug) cout << " -> " << myNoeud->cour ;
				if (autreNoeud != myNoeud->pred)
					throw string (" !ERROR LINKS! ") ;
				if (myNoeud->route != &routes[k][r])
					throw string (" !ERROR ROUTES! ") ;
			}
			if (tracesDebug) cout << endl ;
		}
		if (tracesDebug) cout << endl ;
}

void LocalSearch::reinitSingleDayMoves(Route * r)
{
	for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
		r->nodeAndRouteTested[i] = false ;

	for (Node * tempNoeud = r->depot->nextNode ; !tempNoeud->isDepot ; tempNoeud = tempNoeud->nextNode)
		for (int route = 0 ; route < params->nbVehiclesPerDepot ; route ++)
			depots[tempNoeud->day][route].route->nodeAndRouteTested[tempNoeud->cour] = false ;
}

void LocalSearch::reinitAllSingleDayMoves()
{
	for (int k=1 ; k<= params->nbDays ; k++)
		for (int route = 0 ; route < params->nbVehiclesPerDepot ; route ++)
			for (int i=0 ; i < params->nbClients + params->nbDepots ; i ++ )
				depots[k][route].route->nodeAndRouteTested[i] = false ;
}

bool compPredicateEC(EC_element * i,EC_element * j)
{
	return (i->cost < j->cost - 0.0001 || (i->cost < j->cost + 0.0001 && i->nbCustNodesInChain > j->nbCustNodesInChain));
}

bool LocalSearch::ejectionChains (int day)
{
	SeqData * seq = depots[day][0].seq0_i ;
	int myNodeIndex ;
	int myRouteIndex ;
	Node * myNoeud ;
	Vehicle * myVehicle ;
	Node * myDepot ;
	double myPrevCost ;
	double myDeltaCost ;

	// 1) CHOOSE AN ORDER FOR THE BINS
	vector < Node * > ordreBins ;
	for (int route = 0 ; route < params->nbVehiclesPerDepot ; route ++)
		ordreBins.push_back(&depots[day][route]);
	std::random_shuffle(ordreBins.begin(),ordreBins.end());

	// 2) UPDATE THE DATA STRUCTURE WILL ALL NECESSARY INFORMATIONS AND THE GOOD SIZE
	myRouteIndex = 0 ;
	for (int route = 0 ; route < params->nbVehiclesPerDepot ; route ++)
	{
		myNoeud = ordreBins[route] ;
		// if the route is not empty, this makes a new layer in the ejection chains graph
		if (!myNoeud->nextNode->isDepot)
		{
			ejectionGraph[myRouteIndex][0].myNode = myNoeud ;
			ejectionGraph[myRouteIndex][0].cost = 1.e20 ;
			ejectionGraph[myRouteIndex][0].pred = NULL ;
			ejectionGraph[myRouteIndex][0].nbCustNodesInChain = 0 ;
			ejectionGraph[myRouteIndex][0].routeID = route ;
			myNoeud = myNoeud->nextNode ;
			myNodeIndex = 1 ;
			while (!myNoeud->isDepot)
			{
				ejectionGraph[myRouteIndex][myNodeIndex].myNode = myNoeud ;
				ejectionGraph[myRouteIndex][myNodeIndex].cost = 1.e20 ;
				ejectionGraph[myRouteIndex][myNodeIndex].pred = NULL ;
				ejectionGraph[myRouteIndex][myNodeIndex].nbCustNodesInChain = 0 ;
				ejectionGraph[myRouteIndex][myNodeIndex].routeID = route ;
				myNoeud = myNoeud->nextNode ;
				myNodeIndex ++ ;
			}
			ec_nbElements[myRouteIndex] = myNodeIndex ;
			myRouteIndex ++ ;
		}
	}
	ec_nbRoutes = myRouteIndex ;

	if (ec_nbRoutes < 2)
	{
		//cout << "Only one route for EC, abort" << endl ;
		return false ;
	}

	// 3) SOLVE THE SHORTEST PATH PROBLEM
	// for each bin in order
	for (int r=0 ; r < ec_nbRoutes ; r++)
	{
		myDepot = ejectionGraph[r][0].myNode ;
		myVehicle = myDepot->route->vehicle ;
		myPrevCost = myDepot->pred->seq0_i->evaluation(myDepot->pred->seq0_i,myVehicle); 

		for (int rPrev=0 ; rPrev < r  ; rPrev ++)
		{
			// For the depot node, propagate from another previous depot node
			if (ejectionGraph[rPrev][0].cost < ejectionGraph[r][0].cost - 0.0001 || 
				(ejectionGraph[rPrev][0].cost < ejectionGraph[r][0].cost + 0.0001 && ejectionGraph[rPrev][0].nbCustNodesInChain > ejectionGraph[r][0].nbCustNodesInChain))
			{
				ejectionGraph[r][0].cost = ejectionGraph[rPrev][0].cost ;
				ejectionGraph[r][0].pred = &ejectionGraph[rPrev][0] ;
				ejectionGraph[r][0].bestInsertionPlace = NULL ;
				ejectionGraph[r][0].nbCustNodesInChain = ejectionGraph[rPrev][0].nbCustNodesInChain ;
			}

			// for the depot node, propagate from any predecessor node
			for (int iiPrev=1 ; iiPrev < ec_nbElements[rPrev] ; iiPrev++)
			{
				// for any such predecessor, test the best place of insertion in the route
				for (int ii=0 ; ii < ec_nbElements[r] ; ii++)
				{
					myDeltaCost = seq->evaluation(ejectionGraph[r][ii].myNode->seq0_i,
						ejectionGraph[rPrev][iiPrev].myNode->seq1,
						ejectionGraph[r][ii].myNode->nextNode->seqi_n,
						myVehicle) - myPrevCost ; 
					if (myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][0].cost - 0.0001 || 
						(myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][0].cost + 0.0001 && ejectionGraph[rPrev][iiPrev].nbCustNodesInChain > ejectionGraph[r][0].nbCustNodesInChain))
					{
						ejectionGraph[r][0].cost = myDeltaCost + ejectionGraph[rPrev][iiPrev].cost ;
						ejectionGraph[r][0].pred = &ejectionGraph[rPrev][iiPrev] ;
						ejectionGraph[r][0].bestInsertionPlace = ejectionGraph[r][ii].myNode ;
						ejectionGraph[r][0].nbCustNodesInChain = ejectionGraph[rPrev][iiPrev].nbCustNodesInChain ;
					}
				}
			}
		}

		for (int ii=1 ; ii < ec_nbElements[r] ; ii++)
		{
			// for each customer node in the bin, propagate from the general 0 node
			myDeltaCost = seq->evaluation(ejectionGraph[r][ii].myNode->pred->seq0_i,
				ejectionGraph[r][ii].myNode->nextNode->seqi_n,
				myVehicle) - myPrevCost ;
			if (myDeltaCost < ejectionGraph[r][ii].cost - 0.0001)
			{
				ejectionGraph[r][ii].cost = myDeltaCost ;
				ejectionGraph[r][ii].pred = NULL ;
				ejectionGraph[r][ii].bestInsertionPlace = NULL ;
				ejectionGraph[r][ii].nbCustNodesInChain = 1 ;
			}

			// for each customer node in the bin, propagate from any previous 0 (depot) node.
			for (int rPrev=0 ; rPrev < r  ; rPrev ++)
			{
				if (myDeltaCost + ejectionGraph[rPrev][0].cost < ejectionGraph[r][ii].cost - 0.0001 ||
					(myDeltaCost + ejectionGraph[rPrev][0].cost < ejectionGraph[r][ii].cost + 0.0001 && ejectionGraph[rPrev][0].nbCustNodesInChain >= ejectionGraph[r][ii].nbCustNodesInChain))
				{
					ejectionGraph[r][ii].cost = myDeltaCost + ejectionGraph[rPrev][0].cost ;
					ejectionGraph[r][ii].pred = &ejectionGraph[rPrev][0] ;
					ejectionGraph[r][ii].bestInsertionPlace = NULL ;
					ejectionGraph[r][ii].nbCustNodesInChain = ejectionGraph[rPrev][0].nbCustNodesInChain + 1 ;
				}
			}

			// for each customer node in the bin, propagate from a previous customer node.
			for (int rPrev=0 ; rPrev < r  ; rPrev ++)
			{
				for (int iiPrev=1 ; iiPrev < ec_nbElements[rPrev] ; iiPrev++)
				{
					myDeltaCost = seq->evaluation(ejectionGraph[r][ii].myNode->pred->seq0_i,
						ejectionGraph[rPrev][iiPrev].myNode->seq1,
						ejectionGraph[r][ii].myNode->nextNode->seqi_n,
						myVehicle) - myPrevCost ;

					if (myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][ii].cost - 0.0001 ||
						(myDeltaCost + ejectionGraph[rPrev][iiPrev].cost < ejectionGraph[r][ii].cost + 0.0001 &&  ejectionGraph[rPrev][iiPrev].nbCustNodesInChain >= ejectionGraph[r][ii].nbCustNodesInChain))
					{
						ejectionGraph[r][ii].cost = myDeltaCost + ejectionGraph[rPrev][iiPrev].cost ;
						ejectionGraph[r][ii].pred = &ejectionGraph[rPrev][iiPrev] ;
						ejectionGraph[r][ii].bestInsertionPlace = NULL ;
						ejectionGraph[r][ii].nbCustNodesInChain = ejectionGraph[rPrev][iiPrev].nbCustNodesInChain + 1 ;
					}
				}
			}
		}
	}

	// 4) KEEP THE BEST SHORTEST PATHS in increasing order of costs.
	vector <EC_element *> orderEnds ;
	for (int r=1 ; r < ec_nbRoutes ; r++)
		orderEnds.push_back(&ejectionGraph[r][0]);
	std::sort(orderEnds.begin(),orderEnds.end(),compPredicateEC);

	// 5) APPLY THE MOVES TO UPDATE THE SOLUTION (track back the nodes and apply successive relocates).
	// Apply the best chain
	if (orderEnds[0]->cost < -0.0001)
	{
		EC_element * elementCour ;
		EC_element * elementPred ;
		cout << "----------------- EJECTION CHAINS PHASE ---------------" << endl ;
		for (int myID = 0 ; myID < 1 ; myID++)
		{
			nbEjectionChains ++ ;
			nbEjectionChainsNodes += orderEnds[myID]->nbCustNodesInChain ;
			elementCour = orderEnds[myID];
			elementPred = elementCour->pred ;
			Node * insertionPosition ;
			Node * insertionPositionTemp ;
			while (elementCour != NULL)
			{
				if (!(elementPred == NULL || elementPred->myNode->isDepot))
				{
					//before is a node, after is a depot (insert in the best location, and keep track of where he was.)
					if (elementCour->myNode->isDepot)
					{
						insertionPosition = elementPred->myNode->pred ;		
						insertNoeud(elementPred->myNode,elementCour->bestInsertionPlace);
						reinitSingleDayMoves(elementCour->bestInsertionPlace->route);
						reinitSingleDayMoves(insertionPosition->route);
						elementCour->bestInsertionPlace->route->updateRouteData(false);
						insertionPosition->route->updateRouteData(false);
					}
					else // before is a node, after is a node
					{
						insertionPositionTemp = elementPred->myNode->pred ;
						insertNoeud(elementPred->myNode,insertionPosition);
						reinitSingleDayMoves(insertionPosition->route);
						reinitSingleDayMoves(insertionPositionTemp->route);
						insertionPosition->route->updateRouteData(false);
						insertionPositionTemp->route->updateRouteData(false);
						insertionPosition = insertionPositionTemp ;
					}
				}
				elementCour = elementPred ;
				if (elementPred != NULL) elementPred = elementCour->pred ;
			}
		}
        setRouteEmpty(day);
		return true ;
	}
	else
		return false ;
}

void LocalSearch::addSeqDataInPieces (Node * node, int length, int day)
{
	Node * cour = node ;
	SeqData * courSeq ;
	int courLenght = length + 1 ;
	int temp ;
	while (courLenght > 0)
	{
		temp = min(courLenght,params->sizeSD) ;
		courSeq = cour->seqi_j[temp-1];
		mySeqs.push_back(courSeq);
		courLenght -= temp ;
		if (courLenght != 0) cour = clients[day][courSeq->lastNode].nextNode ;
	}
}

void LocalSearch::addReverseSeqDataInPieces (Node * node, int length, int day)
{
	Node * cour = node ;
	SeqData * courSeq ;
	list <SeqData*> myseqTemp ;
	int courLenght = length + 1 ;
	int temp ;
	while (courLenght > 0)
	{
		temp = min(courLenght,params->sizeSD) ;
		courSeq = cour->seqj_i[temp-1];
		myseqTemp.push_back(courSeq);
		courLenght -= temp ;
		if (courLenght != 0) cour = clients[day][courSeq->firstNode].nextNode ;
	}

	for ( list<SeqData*>::reverse_iterator it=myseqTemp.rbegin() ; it != myseqTemp.rend() ; it++)
		mySeqs.push_back(*it);
}

LocalSearch::LocalSearch(void)
{
	allAttributesSet = false ;
    seqStart = NULL ;
}

// 预处理结构汇总在一个连续区域，但是效益并不高，复现可以考虑不连续的持有
LocalSearch::LocalSearch(Params * params, Individual * individu) : params (params), individu(individu)
{
	nbDays = params->nbDays ;
	int nbVeh ;
	allAttributesSet = true ;
	nbInterSwap = 0 ;
	nbIntraSwap = 0 ;
	nbInter2Opt = 0 ;
	nbIntra2Opt = 0 ;
	nbEjectionChains = 0 ;
	nbEjectionChainsNodes = 0 ;
    seqStart = NULL ;
	nbTotalRISinceBeginning = 0 ;
	nbTotalPISinceBeginning = 0 ;
	vector < Node * > tempNoeud ;
	vector <int> temp2 ;

	clients = new Node * [params->nbDays + 1] ;
	depots = new Node * [params->nbDays + 1] ;
	depotsFin = new Node * [params->nbDays + 1] ;
	routes = new Route * [params->nbDays+1] ;

	for (int kk = 1 ; kk <= params->nbDays ; kk++)
	{
		nbVeh = params->numberVehicle[kk] ;
		clients[kk] = new Node [params->nbClients + params->nbDepots + 1] ;
		depots[kk] = new Node [nbVeh] ;
		depotsFin[kk] = new Node [nbVeh] ;
		routes[kk] = new Route [nbVeh] ; 

		for (int i = 0 ; i <  params->nbClients + params->nbDepots ; i ++ )
			clients[kk][i] = Node(false, i, kk, false, NULL, NULL, NULL, params);

		for (int i = 0 ; i < nbVeh ; i++ )
		{
			depots[kk][i] = Node(true, params->orderVehicles[kk][i].depotNumber, kk, false, NULL, NULL, NULL, params);
			depotsFin[kk][i] = Node(true, params->orderVehicles[kk][i].depotNumber, kk, false, NULL, NULL, NULL, params);
			routes[kk][i] = Route(i, 0, &params->orderVehicles[kk][i], params, individu, kk);
			depots[kk][i].route = &routes[kk][i] ;
			depotsFin[kk][i].route = &routes[kk][i] ;
			routes[kk][i].depot = &depots[kk][i] ;
		}
	}

	for (int i=0 ; i < 4 ; i++)
	{
		resultMoves.push_back(vector<double>(4));
		shouldBeTested.push_back(vector<bool>(4));
	}

	for (int day = 0 ; day <= params->nbDays ; day++) {
        routeOrder.push_back(temp2);
        routeEmpty.push_back(NULL);
    }

    for (int i = params->nbDepots; i < params->nbDepots + params->nbClients; i++)
        routeOrder[0].push_back(i);

	// Initialization for ejection chains
	for (int v=0 ; v < params->nbVehiclesPerDepot ; v ++)
		ejectionGraph.push_back(vector <EC_element> (params->nbClients+1));
	ec_nbElements = vector <int> (params->nbVehiclesPerDepot);

	// Initialization of the SeqData structures, for all nodes
	// These structures were designed to stay contiguous in memory (however it had only a negligible impact on performance)
	// seqStart包含了一个大的seqData数组，所有client的预处理子串连续存储在这个区域
	// 这部分代码用每个node对应的client的index初始化了所有其持有的seqData，并指向对应的数组
	int nbSeqsSet = 0 ;
	int taillemyseqDatas = (params->sizeSD + params->sizeSD + 4) * (params->nbClients + params->nbDepots + 2*params->nbVehiclesPerDepot) * params->nbDays ;
	SeqData * myseqDatas = new SeqData [taillemyseqDatas] ;
    seqStart = myseqDatas ;
	for (int k=1 ; k <= params->nbDays ; k++)
	{
		for (int i=0 ; i < params->nbClients + params->nbDepots ; i++) 
		{
			for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
				myseqDatas[ii].initialisation(clients[k][i].cour,params,individu,k,false);  

			clients[k][i].seq0_i = &myseqDatas[nbSeqsSet] ;
			clients[k][i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
			clients[k][i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
			clients[k][i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				clients[k][i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				clients[k][i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
			clients[k][i].setRemaining();
			nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
		}
        // 分配仓点node的预处理seqData
		for (int i=0 ; i < params->nbVehiclesPerDepot  ; i++)
		{
			for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
				myseqDatas[ii].initialisation(depots[k][i].cour,params,individu,k,false);

			depots[k][i].seq0_i = &myseqDatas[nbSeqsSet] ;
			depots[k][i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
			depots[k][i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
			depots[k][i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depots[k][i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depots[k][i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
			depots[k][i].setRemaining();
			nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
		}
        // depotFin和depot不共用预处理子串，上下两段代码的唯一区别在于depots和depotsFin
		for (int i=0 ; i < params->nbVehiclesPerDepot ; i++)
		{
			for (int ii=nbSeqsSet ; ii < nbSeqsSet+4+params->sizeSD+params->sizeSD ; ii++)
				myseqDatas[ii].initialisation(depotsFin[k][i].cour,params,individu,k,false);

			depotsFin[k][i].seq0_i = &myseqDatas[nbSeqsSet] ;
			depotsFin[k][i].seqi_n = &myseqDatas[nbSeqsSet+1] ;
			depotsFin[k][i].seqi_0 = &myseqDatas[nbSeqsSet+2] ;
			depotsFin[k][i].seqn_i = &myseqDatas[nbSeqsSet+3] ;
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depotsFin[k][i].seqi_j.push_back(&myseqDatas[nbSeqsSet+4+ii]);
			for (int ii=0 ; ii<params->sizeSD ; ii++)
				depotsFin[k][i].seqj_i.push_back(&myseqDatas[nbSeqsSet+4+params->sizeSD+ii]);
			depotsFin[k][i].setRemaining();
			nbSeqsSet += 4+params->sizeSD + params->sizeSD ;
		}
	}
}

LocalSearch::~LocalSearch(void)
{	
	SeqData * seqdeb2 = (SeqData*) seqStart ;
	delete [] seqdeb2 ;

	if (allAttributesSet)
	{
		for (int kk = 1 ; kk <= nbDays ; kk++)
		{
			delete [] clients[kk] ;
			delete [] depots[kk] ;
			delete [] depotsFin[kk] ;
			delete [] routes[kk] ;
		}
		delete [] clients ;
		delete [] depots ;
		delete [] depotsFin ;
		delete [] routes ;
	}
}

