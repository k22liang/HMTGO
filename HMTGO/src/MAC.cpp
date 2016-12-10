// [[Rcpp::plugins(cpp11)]]
#include <stdlib.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
using namespace Rcpp;


NumericVector getPath2TopC(int id, NumericVector parentIndex){
	int pa;
  std::vector<int> path;
	pa=parentIndex[id];
	path.push_back(pa);
	while(pa!=1)
	{
		pa=parentIndex[pa-1];
		path.push_back(pa);
	}
	return(wrap(path));
}

bool not_include(int i, NumericVector path)
{
	for(int j=0;j<path.size();j++)
	{
		if(i==path[j])
		{
			return(false);
		}
	}
	return(true);
}

int maxdeepnode(NumericVector x, NumericVector depth)
{
  int temp=0;
  for(int i=1; i<x.size();i++)
  {
    if(depth[x[i]-1]>depth[x[temp]-1])
    {
      temp=i;
    }
    else if(depth[x[i]-1]==depth[x[temp]-1])
    {
      if (x[i]>x[temp])
      {
        temp=i;
      }
    }
  }
  return(x[temp]);
}

// [[Rcpp::export]]
NumericMatrix MAC(List treedata)
{
	int nTree=as<int>(treedata["nTree"]);
	NumericVector parentIndex = as<NumericVector>(treedata["parentIndex"]);
	NumericVector depth = as<NumericVector>(treedata["depth"]);
	NumericMatrix MA(nTree,nTree);
	std::fill(MA.begin(), MA.end(), -1);
	for(int i=1; i<nTree-1; i++)
	{
	  NumericVector path1=getPath2TopC(i,parentIndex);
		for(int j=i+1; j<nTree;j++)
		{
		  NumericVector path2=getPath2TopC(j,parentIndex);
			if(not_include(i+1,path2))
			{
				NumericVector temp=intersect(path1,path2);
				MA(i,j)=maxdeepnode(temp,depth);
			}
				
		}
	}
	return(MA);
}
