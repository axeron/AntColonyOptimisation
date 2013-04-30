package uk.ac.aber.dcs.ans23.AntColonyOptimisation;

import java.util.*;

/**
 * Class to hold the implementation of an ant colony.
 * @author Andreas O'Brien Svingeseth (ans23@aber.ac.uk)
 * @version 0.1
 */
public class AntColony {
	  protected MinimumSpanningTree mst; //getting the minimum spanning tree from its class.
	  protected double[][] distanceBetweenVertices; //holds the distances between the vertices.
	  protected double[][] nearnessFactor; //holds the nearness factor.
	  protected double[][] pheromoneOnTrail; //holds the pheromone values left on the edges.
	  protected double[][] pheromoneDelta; //hold the delta in the pheromone values.
	  protected double[][] qualityOfEdges; //holds the quality of the edges
	  protected boolean[]  isVisited; //holds what vertices have been visited
	  protected double currentTourLength; //holds the length of the tour of the current ant.
	  protected double bestRunLength; //holds the length of the best tour in the current run.
	  protected double bestTourLength; //holds the length of the best tour so far.
	  
	  protected int[] currentTour; //holds the tour of the current ant 
	  protected int[] bestRun; //holds the best tour in the current run
	  protected int[] bestTour; //holds the best tour so far.
	  
	  private int[] destinations; //holds a buffer for the destinations
	  private double[] sums; //holds a buffer for the probabilities.
	  
	  protected double averageLength; //holds the average tour length.
	  protected double exploitProb; //probability of exploiting the best edge.
	  protected double alpha; //exponent for trail values.
	  protected double beta; //exponent for distances.
	  protected double evaporatinFactor; //holds the pheromone evaporation factor.
	  protected double layexp; //holds the pheromone laying exponent.
	  protected double elite; //holds the best tour enhancement factor.
	  protected int numberOfAnts; //number of ants in current run
	  protected int epochs; //current number of epochs
	  
	  private double maximumTrailValue; //holds the maximum trail value. 
	  private double averageTrailValue; //holds the average trail value. 
	  private Random rand; //random number generator.

	  /**
	   * Constructor. Allows the creation of a new ant colony.
	   * @param mst
	   * @param numberOfAnts
	   * @param rand
	   */
	  public AntColony (MinimumSpanningTree mst, int numberOfAnts, Random rand){ 
	    this.mst = mst; 
	    this.distanceBetweenVertices = mst.distanceBetweenVertices;   
	    int size = mst.size(); //create the trail matrices.
	    this.nearnessFactor = new double[size][size];
	    this.pheromoneOnTrail = new double[size][size];
	    this.pheromoneOnTrail = new double[size][size];
	    this.qualityOfEdges = new double[size][size];
	    this.isVisited = new boolean[size];
	    this.currentTour = new int[size]; 
	    this.currentTourLength = Double.MAX_VALUE;
	    this.bestRun    = new int[size]; 
	    this.bestRunLength   = Double.MAX_VALUE;
	    this.bestTour    = new int[size]; 
	    this.bestTourLength = Double.MAX_VALUE;
	    this.destinations= new int[size];
	    this.sums = new double[size];
	    this.rand = rand; //store the random number generator.
	    this.numberOfAnts  = numberOfAnts; //store the number of ants.
	    this.exploitProb = 0.0; //initializing the probability for exploiting best edge.
	    this.alpha = 1.0; //initializing the weighting exponent for trail values.
	    this.beta = 1.0; //initializing the weighting exponent for distance values.
	    this.evaporatinFactor = 0.1; //initializing the evaporation factor.
	    this.epochs = 0; //initializing the epoch counter.
	  } 

	  /**
	   * Constructor. Allows for the construction of a new ant colony.
	   * @param mst
	   * @param numberOfAnts
	   */
	  public AntColony (MinimumSpanningTree mst, int numberOfAnts){ 
		  this(mst, numberOfAnts, new Random()); 
	  }

	  /**
	   * Allows the program to get the minimum spanning tree problem.
	   * @return the mst problem.
	   */
	  public MinimumSpanningTree getMST (){ 
		  return this.mst; 
	  }

	  /**
	   * Allows the program to set the exploitation probability.
	   * @param exploit
	   */
	  public void setExploitationProbability(double exploit){ 
		  this.exploitProb = exploit;   
	  }
	  
	  /**
	   * Allows the program to set the alpha exponent.
	   * @param alpha
	   */
	  public void setAlpha(double alpha){ 
		  this.alpha = alpha;   
	  }
	  
	  /**
	   * Allows the program to set the beta exponent.
	   * @param beta
	   */
	  public void setBeta(double beta){ 
		  this.beta = beta; 
	  }
	  
	  /**
	   * Allows the program to set the evaporation factor.
	   * @param evap
	   */
	  public void setEvapaporationFactor(double evap){ 
		  this.evaporatinFactor = evap;   
	  }
	  
	  /**
	   * Allows the program to set the trail.
	   * @param trail
	   */
	  public void setTrail(double trail){ 
		  this.elite = elite;   
	  }
	  
	  /**
	   * Allows the program to set the tour enhancement factor.
	   * @param elite
	   */
	  public void setElite(double elite){
		  this.elite = elite;   
	  }

	  /**
	   * Allows the program to get the distances between the vertices
	   * @param i
	   * @param j
	   * @return the distances between the vertices.
	   */
	  public double getDistanceBetweenVertices(int i, int j){
		  return this.distanceBetweenVertices[i][j]; 
	  }
	  
	  /**
	   * Allows the program to get the pheromone values on the trail.
	   * @param i
	   * @param j
	   * @return the pheromone values.
	   */
	  public double getPheromoneValues(int i, int j){
		  return this.pheromoneOnTrail[i][j];  
	  }
	  
	  /**
	   * Allows the program to get the average value trail values.
	   * @return the average trail values.
	   */
	  public double getTrailAverage(){ 
		  return this.averageTrailValue; 
	  }
	  
	  /**
	   * Allows the program to get the maximum trail values.
	   * @return the maximum trail values.
	   */
	  public double getTrailMax(){ 
		  return this.maximumTrailValue;
	  }
	  
	  /**
	   * Allows the program to get the current best tour.
	   * @return the best tour.
	   */
	  public int[] getBestTour(){ 
		  return this.bestTour; 
	  }
	  
	  /**
	   * Allows the program to get the current best tours length.
	   * @return the length of the current best tour.
	   */
	  public double getBestTourLength(){ 
		  return this.bestTourLength;
	  }
	  
	  /**
	   * Allows the program to keep count at what epoch it is currently at.
	   * @return the number of epochs.
	   */
	  public int getEpochs(){ 
		  return this.epochs; 
	  }

	  /**
	   * Allows the program to initialize a new run.
	   */
	  public void init (){ 
		  this.init(-1);   
	  }

	  /**
	   * Allows the program to initialize a new run.
	   * @param val
	   */
	  public void init (double val)
	  {                             /* --- initialize nearness and trail */
	    int    i, j;                /* loop variables */
	    double sum = 0;             /* sum of edge lengths */

	    for (i = this.currentTour.length; --i >= 0; ){
	      for (j = this.currentTour.length; --j >= 0; ){
	        sum += this.distanceBetweenVertices[i][j];/* compute the average tour length */
	      }
	    }
	    
	    this.averageLength = sum /this.currentTour.length;
	    
	    if (val <= 0) val = 1;      /* check and adapt initial value */
	    for (i = this.currentTour.length; --i >= 0; ) {
	      for (j = this.currentTour.length; --j >= 0; ) {
	        this.nearnessFactor[i][j] = Math.pow(this.distanceBetweenVertices[i][j], -this.beta);
	        this.pheromoneOnTrail[i][j] = val; /* compute nearness from distance */
	      }                         /* and set all trail elements */
	    }                           /* to the same value */
	    this.maximumTrailValue = this.averageTrailValue = val;  /* set maximal/average trail value */
	    this.bestTourLength = Double.MAX_VALUE;
	    this.epochs  = 0; 
	  } 



	  private static int find (double vec[], int n, double val)
	  {                             /* --- find edge based on random val. */
	    int i, k;                   /* left and middle element */

	    for (--n, i = 0; i < n; ) { /* do a binary search */
	      k = (i +n) >> 1;          /* in the given vector */
	      if (vec[k] < val) i = k+1;
	      else              n = k;  /* i and n are the boundaries */
	    }                           /* of the section still to search */
	    return i;                   /* return index i for which it is */
	  }  /* find() */               /* vec[i-1] < val <= vec[i] */

	  /*------------------------------------------------------------------*/

	  private void placePheromone (int[] tour, double amount)
	  {                             /* --- place pheromone on ant's tour */
	    int     i;                  /* loop variable */
	    int     src, dst;           /* source and destination of an edge */
	    boolean sym;                /* whether TSP is symmetric */

	    src = this.currentTour[0];         /* get the start of the tour */
	    for (i = tour.length; --i >= 0; ) {
	      dst = src; src = tour[i]; /* traverse the vertices on the tour */
	      this.pheromoneDelta[src][dst] += amount;
	    }                           /* place pheromone on the edge */
	  }  /* placeTrail() */

	  /*------------------------------------------------------------------*/

	  public double runAnt ()
	  {                             /* --- run one ant of the colony */
	    int    i, j, n;             /* loop variables, counter */
	    int    src, dst = 0;        /* source and dest. of next edge */
	    double emax;                /* maximal value of an edge */
	    double sum;                 /* sum of edge values */
	    double chg;                 /* change of trail */

	    /* --- initialize variables --- */
	    for (i = this.isVisited.length; --i >= 0; )
	      this.isVisited[i] = false;  /* clear the visited flags */
	    this.currentTourLength = 0;               /* and the tour length */

	    /* --- run the ant --- */
	    src = this.rand.nextInt(this.currentTour.length);
	    this.currentTour[0]      = src;    /* randomly select the vertex */
	    this.isVisited[src] = true;   /* the ant starts from */
	    for (i = 0; ++i < this.currentTour.length; ) {
	      if ((this.exploitProb > 0)    /* if to exploit best known edge */
	      &&  (rand.nextDouble() < this.exploitProb)) {
	        emax = -1.0;            /* init. the best edge value */
	        for (j = this.currentTour.length; --j >= 0; ) {
	          if (!this.isVisited[j]  /* traverse edges to unvisited verts. */
	          &&  (this.qualityOfEdges[src][j] > emax)) {
	            emax = this.qualityOfEdges[src][j]; dst = j; }
		} }                     /* find the best edge to follow */
	      else {                    /* if to choose edge randomly */
	        sum = 0;                /* init. the quality sum */
	        for (j = n = 0; j < this.currentTour.length; j++) {
	          if (this.isVisited[j]) continue;
	          sum += this.qualityOfEdges[src][j];
	          this.sums[n  ] = sum; /* collect and sum the qualities */
	          this.destinations[n++] = j;   /* of the different edges */
	        }                       /* to unvisited destinations */
	        j   = find(this.sums, n, sum *rand.nextDouble());
	        dst = this.destinations[j];     /* choose destination randomly */
	      }                         /* based on the edge qualities */
	      this.isVisited[dst] = true; /* mark the destination as visited */
	      this.currentTourLength += this.distanceBetweenVertices[src][dst];
	      this.currentTour[i] = src = dst; /* sum the edge lengths and */
	    }                           /* add the vertex to the tour */
	    this.currentTourLength += this.distanceBetweenVertices[src][this.currentTour[0]];

	    /* --- place pheromone --- */
	    chg = this.averageLength/this.currentTourLength; /* compute amount of pheromone */
	    if (this.layexp != 1) chg = Math.pow(chg, this.layexp);
	    this.placePheromone(this.currentTour, chg);

	    return this.currentTourLength;            /* return the length of the tour */
	  } 


	  public double runAllAnts ()
	  {                             /* --- run all ants of the colony */
	    int    i, j;                /* loop variables */
	    double t, min;              /* new/minimal trail value */
	    double stick;               /* stick factor for pheromone */

	    /* --- initialize the edge qualities --- */
	    for (i = this.pheromoneDelta.length; --i >= 0; ) {
	      for (j = this.pheromoneDelta.length; --j >= 0; ) {
	        this.pheromoneDelta[i][j] = 0;   /* init. the trail change matrix */
	        this.qualityOfEdges[i][j] = this.nearnessFactor[i][j]
	                         * ((this.alpha == 1.0) ? this.pheromoneOnTrail[i][j]
	                         : Math.pow(this.pheromoneOnTrail[i][j], this.alpha));
	      }                         /* compute the current qualities */
	    }                           /* of the different edges */

	    /* --- run the ants --- */
	    this.bestRunLength = Double.MAX_VALUE;
	    for (i = this.numberOfAnts; --i >= 0; ) {
	      this.runAnt();            /* run an ant on the trail */
	      if (this.currentTourLength >= this.bestRunLength) continue;
	      System.arraycopy(this.currentTour, 0, this.bestRun, 0, this.bestRun.length);
	      this.bestRunLength = this.currentTourLength;    /* if the new tour is better than */
	    }                           /* the currently best, replace it */
	    if (this.bestRunLength < this.bestRunLength) {
	      System.arraycopy(this.bestRun, 0, this.bestRun, 0, this.bestRun.length);
	      this.bestRunLength = this.bestRunLength;/* if the best run tour is better */
	    }                           /* than the best, replace it */
	    if (this.elite > 0) {       /* strengthen best tour */
	      t = this.averageLength/this.bestRunLength;
	      if (this.layexp != 1) t = Math.pow(t, this.layexp);
	      this.placePheromone(this.bestRun, this.elite *this.numberOfAnts *t);
	    }                           /* place pheromone on best tour */

	    /* --- update trail matrix --- */
	    min = this.averageTrailValue/this.currentTour.length;
	    this.maximumTrailValue = this.averageTrailValue = 0;    /* reinit. the max./avg. trail value */
	    stick    = 1 -this.evaporatinFactor;    /* and compute stick factor */
	    if (this.mst.isSymmetric()){/* if symmetric distances */
	      for (i = this.pheromoneOnTrail.length; --i >= 0; ) {
	        for (j = i; --j >= 0; ) {
	          t = stick     * this.pheromoneOnTrail[i][j]
	            + this.evaporatinFactor *(this.pheromoneDelta[i][j] +this.pheromoneDelta[j][i]);
	          if (t < min) t = min; /* compute the new trail value */
	          this.pheromoneOnTrail[i][j] =    /* from both edge directions */
	          this.pheromoneOnTrail[j][i] = t; /* and store it symmetrically */
	          if (t > this.maximumTrailValue) this.maximumTrailValue = t;
	          this.averageTrailValue += t;        /* update the trail matrix, */
	        }                       /* find highest trail value, */
	      }                         /* and sum the trail values */
	      this.averageTrailValue /= 0.5 *this.currentTour.length *this.currentTour.length; }
	    else {                      /* if asymmetric distances */
	      for (i = this.pheromoneOnTrail.length; --i >= 0; ) {
	        for (j = this.pheromoneOnTrail.length; --j >= 0; ) {
	          t = stick     *this.pheromoneOnTrail[i][j]
	            + this.evaporatinFactor *this.pheromoneDelta[i][j];
	          if (t < min) t = min; /* compute the new trail value */
	          this.pheromoneOnTrail[i][j] = t; /* and store it asymmetrically */
	          if (t > this.maximumTrailValue) this.maximumTrailValue = t;
	          this.averageTrailValue += t;        /* update the trail matrix, */
	        }                       /* find highest trail value, */
	      }                         /* and sum the trail values */
	      if (this.currentTour.length > 1) /* compute average trail value */
	        this.averageTrailValue /= this.currentTour.length *(this.currentTour.length -1);
	    }

	    this.epochs++;               /* count the run */
	    return this.bestRunLength;        /* return length of best tour */
	  }
} 	
