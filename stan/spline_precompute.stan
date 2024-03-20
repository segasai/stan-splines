
  // get the vector of spacings between nodes
  
  vector spline_geths(vector nodes)
  {
    int n = size(nodes) -1;
    vector[n] hs;
    for (i in 1:n)
      {
        hs[i] = nodes[i+1]-nodes[i];
      }
    return hs;
  }
  // obtain the vector of spline coefficients given the location           
  // of the nodes and values there
  // We are using natural spline definition
  
  vector spline_getcoeffs(vector nodes, vector vals)
  {
    int n_nodes = size(nodes);
    int n=n_nodes-1;
    vector[n] hi;
    vector[n] bi;
    vector[n-1] vi;
    vector[n-1] ui;
    vector[n_nodes] ret;
    vector[n-1] zs;
    matrix[n-1,n-1] M = rep_matrix(0, n-1, n-1);

    n = n_nodes-1;

    for (i in 1:n)
      {
        hi[i] = nodes[i+1]-nodes[i];
        bi[i] =  1/hi[i]*(vals[i+1]-vals[i]);
      }
    for (i in 2:n)
      {
        vi[i-1] = 2*(hi[i-1]+hi[i]);
        ui[i-1] = 6*(bi[i] - bi[i-1]);
      }
    for (i in 1:n-1)
      {
        M[i,i] = vi[i];
      }
    for (i in 1:n-2)
      {
        M[i+1,i] = hi[i+1];
        M[i,i+1] = hi[i+1];
      }
    //print (M)
    zs = M \ ui ; //mdivide_left_spd(M, ui); 
    ret[1]=0;
    ret[n_nodes] =0;
    ret[2:n_nodes-1]=zs;

    return ret;
  }
  // Evaluate the spline, given nodes, values at the nodes             
  // spline coefficients, locations of evaluation points
  // and integer bin ids of each point
  vector spline_eval(vector nodes,
                     vector vals, vector zs,
                     matrix mat, array[] int pos)
  {
    int n_nodes = size(nodes);
    int n_dat = size(pos);
    vector[n_nodes-1] h;
    vector[n_dat] ret;
    vector[n_nodes-1] A; // not quite A but A/6h                                            
    vector[n_nodes-1] B;
    vector[n_nodes-1] C;
    vector[n_nodes-1] D;
    h = spline_geths(nodes);
    for(i in 1:n_nodes-1)
      {
        A[i] = zs[i+1] / 6 / h[i];
        B[i] = zs[i] / 6 / h[i];
        C[i] = vals[i+1] / h[i] - h[i] * zs[i+1] / 6;
        D[i] = vals[i] / h[i] - h[i] * zs[i] / 6;
      }
    ret = (
           A[pos] .* mat[:,1] +
           B[pos] .* mat[:,2] +
           C[pos] .* mat[:,3] +
           D[pos] .* mat[:,4]
           );
    return ret;
  }
  
  matrix spline_getmat(vector x, vector nodes, array[] int pos)
  // obtain a list of vectors
  // (x-nodes[i])^3, (nodes[i+1]-x)^3, (x-nodes[i]), (nodes[i+1]-x)
  {
    int n_nodes = size(nodes);
    int n_dat = size(pos);
    matrix[n_dat,4] mat;
    for (i in 1:n_dat)
      {
	int curpos = pos[i] ;
	int curpos1 = curpos + 1;
	// Filling the matrix with (x-Left) (Right - x) (x-Left)^3 (Right-x)^3
	mat[i, 3] = (x[i]-nodes[curpos]);
	mat[i, 4] = (nodes[curpos1]-x[i]);
	mat[i, 1] = mat[i,3] .* mat[i,3] .* mat[i,3];
	mat[i, 2] = mat[i,4] .* mat[i,4] .* mat[i,4];
      }
    return mat;
  }
  
  // find in which node interval we should place each point of the vector      
  array[] int spline_findpos(vector nodes, vector x)
  {
    int n_dat = size(x);
    int n_nodes = size(nodes);
    int ret[n_dat];
    for (i in 1:n_dat)
      {
	int success = 0;
       	for (j in 1:n_nodes-1)
          {
            if ((x[i]>=nodes[j]) && (x[i]<nodes[j+1]))
              {
                ret[i] = j;
		success=1;
		break;
              }
          }
	if (success==0)
	  {
	    reject("Point outside knot");
	  }
      }
    return ret;
  }
