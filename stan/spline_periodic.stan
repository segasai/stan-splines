
  // get the vector of spacings between nodes
  vector spline_geths(vector nodes)
  {
    int n = size(nodes);
    vector[n-1] hs;
    for (i in 1:n-1)
      {
        hs[i] = nodes[i+1]-nodes[i];
      }
    return hs;
  }
  
  // obtain the vector of spline coefficients given the location
  // of the nodes and values there
  // We are using natural spline definition           
  vector spline_getcoeffs(vector nodes, vector vals0)
  {
    int n = size(nodes);
    vector[n] vals;
    
    vector[n - 1] hi;
    vector[n - 1] bi;
    vector[n - 1] vi;
    vector[n - 1] ui;
    
    vector[n] ret;
    vector[n - 1] zs;
    matrix[n - 1,n - 1] M = rep_matrix(0, n - 1, n - 1);

    for (i in 1:n-1)
      {
	vals[i] =vals0[i];
      }
    vals[n] = vals0[1];
    for (i in 1:n-1)
      {
        hi[i] = nodes[i+1]-nodes[i];
        bi[i] =  1/hi[i]*(vals[i+1]-vals[i]);
      }
    for (i in 1:n-1)
      {
        vi[i] = 2*(hi[i] + hi[i%(n-1)+1]);
        ui[i] = 6*(bi[i%(n-1)+1] - bi[i]);
      }
    for (i in 1:n-1)
      {
        M[i,i] = hi[i];
        M[i,i%(n-1) + 1] = vi[i];
        M[i,(i+1)%(n-1) +1] = hi[i%(n-1)+1];
      }
    zs = M \ ui ;
    ret[1:n-1]=zs;    
    ret[n] = zs[1]; // periodic condition
    return ret;

  }
  // Evaluate the spline, given nodes, values at the nodes
  // spline coefficients, locations of evaluation points
  // and integer bin ids of each point            
  vector spline_eval(vector nodes,
                     vector vals0, vector zs,
                     vector x, array[] int i)
  {
    int n = size(nodes);
    int n_dat = size(x);
    vector[n] vals;
    vector[n-1] h;
    vector[n_dat] ret;
    array[n_dat] int i1;
    for (ii in 1:n_dat)
      {
        i1[ii] = i[ii] + 1;
      }
    for (ii in 1:n-1)
      {
	vals[ii]= vals0[ii];
      }
    vals[n] = vals0[1];
    h = spline_geths(nodes);

    ret = (
           zs[i1] ./ 6 ./ h[i] .* square(x-nodes[i]) .* (x-nodes[i])+
           zs[i]  ./ 6 ./ h[i] .* square(nodes[i1]-x) .* (nodes[i1]-x)+
           (vals[i1] ./ h[i] - h[i] .* zs[i1] ./ 6) .* (x-nodes[i])+
           (vals[i] ./ h[i] - h[i] .* zs[i] ./ 6) .* (nodes[i1]-x)
           );
    return ret;
  }
  // find in which node interval we should place each point of the vector                   
array[] int spline_findpos(vector nodes, vector x)
  {
    int n_nodes = size(nodes);
    int n_dat = size(x);
    array[n_dat] int ret;
    for (i in 1:n_dat)
      {
	int success = 0;
        for (j in 1:n_nodes-1)
          {
            if ((x[i]>=nodes[j]) && (x[i]<nodes[j+1]))
              {
                ret[i] = j;
		success = 1;
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
