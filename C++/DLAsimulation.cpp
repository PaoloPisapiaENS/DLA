///////////////////////////////////////////////////////////////////
//    We want to define a DLA simulator                          //
//    This is the modified DLA:                                  //
//      we consider the vertex whence it attached to the cluster //
//      not just the boundary of the cluster                     // 
///////////////////////////////////////////////////////////////////


// Include a header file containing all the necessary libraries and function prototypes
#include "/mnt/d/Offline_Documents/University/PhD_Paris/PhD_work/Simulations/Simple_Graph/my_UpgradedHeader.h" 

#define PRINT 0 // 0: don't print anything; 1: print everything exept for the RW steps; 2: print everything; 

int g_notArrivedCounter = 0;
// Counters for the total number of attempts 
int g_totalAttemptsCounter = 0;
int g_successfulAttemptsCounter = 0;


/******************************************************************
//                Initialization Function                        //
******************************************************************/
void DLA( Graph  /*Provide a graph to run the random walk on*/,
                          //int* trail  /*Provide a pointer to an array to store the trail of the RW*/,
                          int*  /*Provide the array containing the desired trail of the Laplacian RW*/,
                          int    /*as well as its size*/,
                          int  nSteps=1e6  /*(OPTIONAL) the number of steps (just to make sure the loop ends)*/);

////////////////////////////////////////////
//                MAIN                    //
////////////////////////////////////////////

int main()
{
  cout << setprecision(4); 

  srand48(time(NULL)); // Initialize the sequence

  // We want to compute a mean of the probabilities, so we need a variable to store the numbers of runs
  int nRun = 1;

  // Variable to store the partial sum of the probabilities to compute their mean
  double partProbSum = 0.;

  // Number of iterations
  int N= 5e6;

  // Typical dimension of the graph
  int nbase=3;

  // Shape of the tiles (4=square)
  int v=4;

  // Absorption weight
  double q = 0;

  // Define the condition, i.e. what Laplacian RW we want
  int condition[] = {1,0,4,5,7};


  ///////////////////////////////////////////////////////////////////

  // Let's print out that this is the Laplacian RW
  cout << "\n DLA" << endl;

  // Let's print out the condition
  cout << "\n We want to obtain the cluster: (" << condition[0] << ", ";
  for (int i = 1; i < sizeof(condition)/sizeof(condition[0]) -1 ; i++)
    cout << condition[i] << ", " ;
  cout << condition[sizeof(condition)/sizeof(condition[0]) -1] << ")" << endl;
    

  // Create the graph
  Graph g=RegularGraph(nbase,v, q);

  // Retreive the graph dimension
  int n=g.graph_size();

  // Compute the mean of the probabilities
  for (int j = 0; j < nRun; j++)
  {
    for(int i=0; i<N; i++)
    {
      DLA(g, condition, sizeof(condition)/sizeof(condition[0]) );

      #if PRINT==1 ||  PRINT==2
      cout << "\n******************************************************" << endl;
      #endif
    }

    partProbSum += (double) g_successfulAttemptsCounter/(g_totalAttemptsCounter);

    g_successfulAttemptsCounter = 0;
    g_totalAttemptsCounter = 0;
  }
  
  // Print out the numeber of correct DLA
  cout << "\n In " << nRun << " runs of " << N << " attempts each, the mean probability to get the DLA cluster (" << condition[0] << ", ";
  for (int i = 1; i < sizeof(condition)/sizeof(condition[0]) -1 ; i++)
    cout << condition[i] << ", " ;
  cout << condition[sizeof(condition)/sizeof(condition[0]) -1] << ") is " << partProbSum/nRun * 100 << "%" << endl;

  return 0;

} // End of main()




////////////////////////////
//  Function definitions  //
////////////////////////////

// We make a slight change to the RandomWalk (RW) function in order to make it stop after a fixed number of steps 'nSteps'
void DLA( Graph g /*Provide a graph to run the random walk on*/,
                          //int* trail  /*Provide a pointer to an array to store the trail of the RW*/,
                          int* condition /*Provide the array containing the desired trail of the Laplacian RW*/,
                          int conditionSize    /*as well as its size*/,
                          int  nSteps  /*(OPTIONAL) the number of steps (just to make sure the loop ends)*/)
{  
  // Set the starting vertex of the RWs as the last vertex of the condition
  int v = condition[conditionSize-1];

  // Define the number of steps left and initialize it
  int nStepsLeft = nSteps;

  // Define the random variable to decide where to move
  double x; 

  // Retrieve the size of the graph
  int n = g.graph_size();

  // Total rate buffer (the total rate information is stored in the n+1 column of the adjacency matrix)
  double totR= g.graph_edge(v,n);

  // Total probability buffer
  double totP = 1.;

  // Buffer for the probability to move towards a specific vertex
  double* moveP = new double[n];
  // And set its index variable
  int k = 0;

  // Start by storing the first vertex in the array 'trail'. It will be useful to then performe the loop erasure. 
  //trail[nSteps-nStepsLeft] = v;
  
  // Define the index that will run on the condition array 
  int j=0;
  
  // Define the growth vertex buffer, where the position of the Laplacian RW is at a given moment, initialized as the first vertex of the condition array
  int tip = condition[j];  

  // Define a local attempt counter
  int attemptsCounter = 0;

  // Create a map where the Laplacian RW will be stored. Start by adding the initial point, which must be the 'tip'
  unordered_map<int, bool> mp;

  mp[tip] = 1;

  while(mp.find(condition[conditionSize-1]) == mp.end()) // Runs until the DLA reaches the source
  {
    // Reset the starting vertex as the last vertex in the condition...
    v = condition[conditionSize-1]; 
    
    // ...and the starting rate...
    totR = g.graph_edge(v,n); 

    #if PRINT==1 ||  PRINT==2
      cout << "\n From the initial vertex " << v << endl ;
    #endif 

    while( totR * nStepsLeft !=0) // The loop stops when we can no longer move from a vertex, i.e. when its not connected to any other vertex (and hence the total rate is zero), 
                                  // OR when it has performed all the steps OR when it hits the trail of the Laplacian RW
    {

      // Resets the total probability to 1 at the beginning of each 'while' loop
      totP = 1.;                
      
      // Retrieve the move probabilities as the edge weight divided by the total weight 
      for ( int i=0; i<n; i++)
        moveP[i] = g.graph_edge(v,i)/totR;  

      // Now it initializes the random variable and assign the running index 'j' to start from the end (the last vertex is the absorption one)
      x = drand48();   // Assign a random variable
      k = n;           // Assign the index to start from the end (the last vertex is the absorption one)

      #if PRINT==2
        cout << "\n (The random number is: " << x << ")" << endl ;
      #endif

      // This loop checks in which probability range one finds itself: when it breaks out, the index k represents the edge one moves to
      while ( x<= totP)               
      {
        k--;
        totP -= moveP[k];
      }
      
      #if PRINT==2
        if(k==(n-1))        // i.e. if it gets absorbed (the absorption vertex is the last one)
          cout << "Absorbed at vertex " << v  << endl << endl; // We get absorbed at that vertex
      
        else
          cout << " We move to vertex " << k << " (probability interval: [" << totP << ", " << totP+moveP[k] << "] )";    // We move to the new vertex k

        if (v==k) 
          cout << "\nA SELF-LOOP!!\n";
      #endif

      // If k is not present in condition[], good, go ahed with the RW
      if (mp.find(k) == mp.end())
      {
        #if PRINT==1 ||  PRINT==2
          cout << "\n Vertex " << k << " not in the trail :)";
        #endif

        // Finally, reduce the number of steps left..
        nStepsLeft--;   

        //...update the starting vertex...
        v=k;

        // ...and the starting rate...
        totR = g.graph_edge(v,n); 
              
        // ... and record the starting vertex in the trail array
        //trail[nSteps-nStepsLeft] = v;
      }
      // If instead k is present, i.e. the RW hit the DLA cluster, report it by braking the loop
      else
        goto jump1;
    }
      

    // If we jump here is because we hit the DLA cluster
    jump1:
      // Store the previous vertex in the map
        mp[v] = 1;

      #if PRINT==1 ||  PRINT==2
        cout << "\n Vertex " << k << " in the cluster";
        cout << "\n So add the preceding vertex to the cluster";
        cout << " (after " << attemptsCounter << " failed attempts)" << endl;
      #endif
      // Then, redo the RW again
  }

  // Finally, check if the trail is the same as the condition: first check their size...
  if (conditionSize == mp.size())
  {
    // ...then check whether all the elements of condition[] are in the Laplacian RW just created
    for(int i=0; i<conditionSize; i++ )
    {
      // If the element 'i' is not found, it means it is not the desired cluster. Break the 'for' and go directly to the total attempt counter increment
      if(mp.find(condition[i]) == mp.end()) 
      {
        #if PRINT==1 ||  PRINT==2
          cout << "\n Wrong cluster :(" << endl;
        #endif

        goto jump2;
      }
    }

    // Otherwise, if all the elemnts are present, count it as a succesful attempt      
    g_successfulAttemptsCounter++;
  }
  #if PRINT==1 ||  PRINT==2
    else
      cout << "\n Wrong cluster :(" << endl;
  #endif
  

  jump2:
    g_totalAttemptsCounter++;
    attemptsCounter++;
  
  #if PRINT == 1 ||  PRINT==2  
    cout << "\n DLA cluster obtained in " << attemptsCounter << " attempts" << endl;
  #endif

  delete[] moveP;

} // End of function DLA