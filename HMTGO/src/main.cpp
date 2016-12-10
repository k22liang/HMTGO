#include <stdlib.h>
#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;
using namespace std;




NumericVector dmixtureC(NumericVector x, double alpha, double beta, double p, bool iflog = 0) {
	NumericVector lk;
	if (p == 0){
		return(dbeta(x, alpha, beta, iflog));
	}
	lk = (1 - p)*dbeta(x, alpha, beta, 0) + p;
	if (iflog){
		return(log(lk));
	}
	else{
		return(lk);
	}

}


NumericVector getOptimAltC(NumericVector x, NumericVector pval,
	NumericVector weight, bool penalty, std::string Method, List control_R_alt, Function optim_R, Function likelihoodfunc){
	List result;
	result = optim_R(Named("par", x), Named("fn", likelihoodfunc), Named("pval", pval), Named("weight", weight), Named("penalty", penalty), Named("method", Method), Named("control", control_R_alt));
	int converge;
	NumericVector par(2), res(3);
	converge = as<int>(result["convergence"]);
	par = as<NumericVector>(result["par"]);
	res[0] = par[0];
	res[1] = par[1];
	res[2] = as<double>(result["value"]);

	if (converge == 0){  //converge
		return (res);
	}
	else if (converge == 1){
		//cout << "getOptimAlt used up iteration!" << endl;
		return(res);
	}
	else if (converge == 10){
		//cout << "degeneracy of the Nelder - Mead simplex" << endl;
		result = optim_R(Named("par", x), Named("fn", likelihoodfunc), Named("pval", pval), Named("weight", weight), Named("penalty", penalty), Named("method", "BFGS"), Named("control", control_R_alt));
		converge = as<int>(result["convergence"]);
		par = as<NumericVector>(result["par"]);
		res[0] = par[0];
		res[1] = par[1];
		res[2] = as<double>(result["value"]);
		if (converge == 0){
			return(res);
		}
		else {
			//cout << "Convergence problem in getOptimAlt!" << endl;
			return(0);
		}
	}
	else{
		//cout << "Convergence value is not readable in getOptimAlt!" << endl;
		return(0);
	}

}


NumericVector getOptimNullC(NumericVector x, NumericVector pval,
	NumericVector weight, std::string Method, List control_R_null, Function optim_R, Function likelihoodfunc){
	List result;
	NumericVector pval_temp = clone(pval);
	double m_c;
	sort(pval_temp.begin(), pval_temp.end());
	int k = 1;
	while (pval_temp[k] == pval_temp[0]){
		k += 1;
	}
	m_c = pval_temp[k];

	result = optim_R(Named("par", x), Named("fn", likelihoodfunc), Named("pval", pval), Named("weight", weight), Named("m_c", m_c), Named("method", Method), Named("control", control_R_null));
	int converge;
	NumericVector par(3), res(4);
	converge = as<int>(result["convergence"]);

	if (converge == 10){  //degeneracy of the Nelder_Mead simplex
		//cout << "degeneracy of the Nelder_Mead simplex" << endl;
		result = optim_R(Named("par", x), Named("fn", likelihoodfunc), Named("pval", pval), Named("weight", weight), Named("method", "BFGS"), Named("control", control_R_null));
		converge = as<int>(result["convergence"]);
	}
	par = as<NumericVector>(result["par"]);
	res[0] = par[0];
	res[1] = par[1];
	res[2] = par[2];
	res[3] = as<double>(result["value"]);
	if (converge == 0){
		return(res);
	}
	else {
		if (converge == 1){
			//cout << "getOptimNull used up iteration!" << endl;
			//cout << res << endl;
			return(res);
		}
		else{
			//cout << "Convergence problem in getOptimNull!" << endl;
			//cout << res << endl;
			return(0);
		}
	}

}

double getCumProbC(int x, int mn, NumericVector cp, NumericVector condProd, NumericVector parentIndex)
{
	int target = x;
	double result = cp[x];
	while (target != mn)
	{
		result *= condProd[target];
		target = parentIndex[target] - 1;
	}
	return (result);
}

// [[Rcpp::export]]
NumericVector getPDEC(NumericVector condProb, NumericVector parentIndex, NumericMatrix componMappingMat, List sepLocation, List formula)
{
	int nTree = componMappingMat.ncol(), nDAG = componMappingMat.nrow();
	NumericVector ccpp(nTree), pde(nDAG);
	ccpp = ccpp + condProb[0];
	for (int i = 1; i < nTree; i++)
	{
		ccpp[i] = ccpp[parentIndex[i] - 1] * condProb[i];
	}
	int n_singleCompon = 0;
	NumericVector  singleCompon(nDAG);


	for (int i = 0; i < nDAG; i++)
	{
		double total = 0;
		for (int j = 0; j < nTree; j++) {
			total += componMappingMat(i, j);
		}
		if (total == 1)
		{
			singleCompon[i] = 1;
			n_singleCompon += 1;
		}
		else
		{
			singleCompon[i] = 0;
		}
	}


	NumericVector componList(n_singleCompon);
	int l = 0;
	for (int i = 0; i < nDAG; i++)
	{
		if (singleCompon[i] == 1)
		{
			int k = 0;
			while (componMappingMat(i, k) == 0){
				k += 1;
			}
			componList[l] = k + 1;
			pde[i] = ccpp[k];
			l += 1;
		}
		else
		{

			NumericVector cp(nTree,1.0);
			int start_pos = 0, mn=0;
			NumericVector temp_sepLocation = as<NumericVector>(sepLocation[i]);



			for (int j = 0; j < temp_sepLocation.size(); j++)
			{
				int last_pos = temp_sepLocation[j] - 1;
				NumericVector temp_formula = as<NumericVector>(formula[i]);
				mn = temp_formula[start_pos] - 1;//index
				double total = 1;
				for (int m = 1; m <= (last_pos - start_pos); m++)
				{
					int temp = temp_formula[(start_pos + m)] - 1;//index
					total *= (1 - getCumProbC(temp, mn, cp, condProb, parentIndex));
				}
				cp[mn] = 1 - total;
				start_pos = last_pos + 1;
			}


			if (mn != 1)
			{
				pde[i] = getCumProbC(mn, 0, cp, condProb, parentIndex)*condProb[0];
			}
			else
			{
				pde[i] = cp[0] * condProb[0];
			}

		}
	}

	return(pde);
}

// [[Rcpp::export]]
List HMTEMC(Function optim_R, Function likelihoodfunc_Alt, Function likelihoodfunc_Null, NumericVector init, NumericVector pval,
	NumericVector depth, NumericVector leaf, NumericVector nonleaf, NumericMatrix childrenIndexMat,
	NumericVector propInParent, NumericVector childrenCount, NumericVector parentIndex, int nDAG, double betaTemp = 1, std::string null_dist = "Mixture",
	bool penalty = false, std::string num_method = "Nelder-Mead", int max_it = 200, int innerit = 1000, double conv_threshold = 1e-8){


	double alpha, alpha_last, beta, beta_last, q, q_last, q0, alpha_null, beta_null, lambda;
  //double alpha_null_last, beta_null_last,lambda_last;
	double sum_P11, sum_P10, g1, g0, temp, prod0, prod1;
	int nTree = pval.size(), i_leaf;

	alpha = init[0];
	beta = init[1];
	q = init[2];
	q0 = init[3];

	alpha_null = init[4];
	beta_null = init[5];
	lambda = init[6];

	NumericVector rho(nTree,-1.0), LL(max_it), Beta(nTree,-1.0), BetaTran0(nTree,-1.0), transProb(nTree), P1(nTree), zerovector(nTree);
	NumericVector BetaTran1(nTree,-1.0), Nu(nTree,-1.0), propInParentTemp(nTree), den(nTree), den_null(nTree,1.0);
	LogicalVector iQ(nTree);


	propInParent[0] = 0;
	int counter_record = 0;
	//log likelihood
	for (int counter = 0; counter < max_it; counter++){
		counter_record = counter;
		q_last = q;
		q = pow(q, betaTemp) / (pow(q, betaTemp) + pow(1 - q, betaTemp));
		q0 = pow(q0, betaTemp) / (pow(q0, betaTemp) + pow(1 - q0, betaTemp));
		for (int i = 0; i < nTree; i++){
			propInParentTemp[i] = pow(propInParent[i], betaTemp) / (pow(propInParent[i], betaTemp) + pow(1 - propInParent[i], betaTemp));
			iQ[i] = (q >= propInParentTemp[i]);
		}
		transProb = pmax(q, propInParentTemp);
		transProb[0] = q0;
		P1 = zerovector + q0;
		for (int i = 1; i < nTree; i++){
			P1[i] = P1[(parentIndex[i] - 1)] * transProb[i];
		}


		den = pow(dbeta(pval, alpha, beta), betaTemp);


		if (null_dist == "Beta"){
			den_null = pow(dbeta(pval, alpha_null, beta_null), betaTemp);
		}
		else if (null_dist == "Mixture"){
			den_null = pow(dmixtureC(pval, alpha_null, beta_null, lambda), betaTemp);
		}
		else{ ; }

		//Upward recursion

		for (int i = 0; i < leaf.size(); i++){
			i_leaf = leaf[i] - 1;
			g1 = den[i_leaf] * P1[i_leaf];
			g0 = den_null[i_leaf] * (1 - P1[i_leaf]);
			Nu[i_leaf] = g1 + g0;
			Beta[i_leaf] = g1 / (g0 + g1);
			BetaTran0[i_leaf] = (1 - Beta[i_leaf]) / (1 - P1[i_leaf]);

			if (Beta[i_leaf] == 0){
				temp = 0;
			}
			else{
				temp = Beta[i_leaf] * transProb[i_leaf] / P1[i_leaf];
			}
			BetaTran1[i_leaf] = (1 - Beta[i_leaf])*(1 - transProb[i_leaf]) / (1 - P1[i_leaf]) + temp;
		}


		//nonleaf was built in 2.tree.transformation.R in a bottom - up fashion, thus, we can assume all children have been updated if we follow the order of nonleaf

		for (int i = 0; i < nonleaf.size(); i++){
			i_leaf = nonleaf[i] - 1;
			prod0 = 1;
			prod1 = 1;
			for (int j = 0; j < childrenCount[i_leaf]; j++){
				prod0 *= BetaTran0[childrenIndexMat(j, i_leaf) - 1];
				prod1 *= BetaTran1[childrenIndexMat(j, i_leaf) - 1];
			}
			g0 = prod0*den_null[i_leaf] * (1 - P1[i_leaf]);
			g1 = prod1*den[i_leaf] * P1[i_leaf];
			Nu[i_leaf] = g1 + g0;
			Beta[i_leaf] = g1 / (g0 + g1);
			BetaTran0[i_leaf] = (1 - Beta[i_leaf]) / (1 - P1[i_leaf]);
			if (Beta[i_leaf] == 0){
				temp = 0;
			}
			else{
				temp = Beta[i_leaf] * transProb[i_leaf] / P1[i_leaf];
			}
			BetaTran1[i_leaf] = (1 - Beta[i_leaf])*(1 - transProb[i_leaf]) / (1 - P1[i_leaf]) + temp;

		}

		// the larger the better
		LL[counter] = sum(log(Nu));

		//Downward recursion
		NumericVector P11(nTree), P10(nTree);

		if (q == 0){
			rho = clone(zerovector);
			rho[0] = Beta[0];
			P11 = clone(zerovector);
			for(int h=0; h<nTree; h++)
			{
				P10[h] = zerovector[h] + 1;
			}
			P10[0] = 0;
		}
		else{
			rho[0] = Beta[0];
			for (int i = 1; i < nTree; i++){
				rho[i] = Beta[i] / P1[i] * (rho[parentIndex[i] - 1] * transProb[i] / BetaTran1[i]);
			}
			P11[0] = 0;
			P10[0] = 0;
			for (int i = 1; i < nTree; i++){
				// 2-slice probability in the M step
				// P11 = P(S_i = 1, S_rho(i) = 1 | \p)
				P11[i] = Beta[i] * transProb[i] * rho[parentIndex[i] - 1] / P1[i] / BetaTran1[i];
				//P10 = P(S_i = 0, S_rho(i) = 1 | \p)
				P10[i] = (1 - Beta[i])*(1 - transProb[i]) * rho[parentIndex[i] - 1] / (1 - P1[i]) / BetaTran1[i];
			}
		}

		//use numerical routines to get MLE of beta parameters
		NumericVector result1(3), result2(4), par1(2), par2(3);
		List control_R = List::create(Named("maxit") = innerit);
		par1[0] = alpha;
		par1[1] = beta;
		alpha_last = alpha;
		beta_last = beta;
		result1 = getOptimAltC(par1, pval, rho, penalty, num_method, control_R, optim_R, likelihoodfunc_Alt);

		alpha = result1[0];
		beta = result1[1];

		if (null_dist == "Mixture"){
			par2[0] = alpha_null;
			par2[1] = beta_null;
			par2[2] = lambda;
      //alpha_null_last = alpha_null;
			//beta_null_last = beta_null;
			//lambda_last = lambda;
			result2 = getOptimNullC(par2, pval, 1 - rho, num_method, control_R, optim_R, likelihoodfunc_Null);
			alpha_null = result2[0];
			beta_null = result2[1];
			lambda = result2[2];
		}

		sum_P11 = 0;
		sum_P10 = 0;
		for (int i = 0; i<nTree; i++){
			if (iQ[i] == true){
				sum_P11 += P11[i];
				sum_P10 += P10[i];
			}
		}


		if (sum_P11 == 0){
			q = 0;
		}
		else{
			q = sum_P11 / (sum_P11 + sum_P10);
		}
		q0 = rho[0];
		//return(List::create(Named("test") = q, Named("test1") = iQ, Named("test2") = sum_P11));
		//break;
		if (pow(pow((q - q_last), 2) + pow((alpha - alpha_last), 2) + pow((beta - beta_last), 2), 0.5) < conv_threshold)
		{
			break;
		}
	}//end counter
	NumericVector cp(nTree);
	if (betaTemp == 1){
		if (q == 0 || q0 == 0){
			;
		}
		else{
			cp[0] = rho[0];
			for (int i = 1; i < nTree; i++){
				if (rho[parentIndex[i] - 1] == 0){
					cp[i] = 0;
				}
				else{
					cp[i] = rho[i] / rho[parentIndex[i] - 1];
				}
			}
		}
	}
	if (null_dist == "Mixture"){
		NumericVector par_ret(7);
		par_ret[0] = alpha;
		par_ret[1] = beta;
		par_ret[2] = q;
		par_ret[3] = q0;
		par_ret[4] = alpha_null;
		par_ret[5] = beta_null;
		par_ret[6] = lambda;
		if (betaTemp == 1){
			return List::create(Named("par") = par_ret, Named("LL") = LL[counter_record], Named("iter") = counter_record + 1, Named("cp") = cp);
		}
		else{
			return List::create(Named("par") = par_ret, Named("LL") = LL[counter_record], Named("iter") = counter_record + 1);
		}
	}
	else{
		NumericVector par_ret(4);
		par_ret[0] = alpha;
		par_ret[1] = beta;
		par_ret[2] = q;
		par_ret[3] = q0;
		if (betaTemp == 1){
			return List::create(Named("par") = par_ret, Named("LL") = LL[counter_record], Named("iter") = counter_record + 1, Named("cp") = cp);
		}
		else{
			return List::create(Named("par") = par_ret, Named("LL") = LL[counter_record], Named("iter") = counter_record + 1);
		}
	}

}

// [[Rcpp::export]]
double getNegNullLlikeC(NumericVector x, NumericVector pval, NumericVector weight, double m_c=0){
  if ((x[0] <= 1) || (x[0] > 20) || (x[1] <= 1) || x[1]>20 || (((x[0]-1)/(x[0]+x[1]-2)) < m_c) || (x[2]<0) || (x[2]>1) ){
    return(1e+10);
  }
  return(-sum(weight*dmixtureC(pval, x[0], x[1], x[2], 1)));

}

// [[Rcpp::export]]
double getNegAltLlikeC(NumericVector x, NumericVector pval, NumericVector weight, bool penalty){
  double p_betavalue;
  double temp, penal = 0;
  int nTree = pval.size();
  p_betavalue = R::pbeta(0.5, x[0], x[1], 0, 0);
  if ((x[0] <= 0) || (x[0] > 1) || (x[1] <= 1)){
    return(1e+20);
  }
  if (penalty){
    if (p_betavalue< 0.95){
      penal = (0.95 - p_betavalue)*100*nTree;
    }
    if (p_betavalue < 0.5){
      penal += (0.5 - p_betavalue)*100*nTree;
    }
  }
  temp = sum(weight*dbeta(pval, x[0], x[1], 1));
  penal = penal- temp;
  return(penal);
}

// [[Rcpp::export]]
List HMT_DAEMC(NumericVector pval, List treedata, Function optim_R, Function likelihoodfunc_Alt, Function likelihoodfunc_Null, NumericVector schedule,
	NumericVector init_origin, bool penalty = false, std::string null_dist = "Mixture", std::string num_method = "Nelder-Mead", int iternum_final = 1000,
	int iternum_hot = 100, int iternum_inner=1000, double conv_threshold = 1e-8)
{
	double betaTemp;
	List result;



	std::string logfile = "EM_log.txt";
	NumericVector depth, leaf, nonleaf, childrenCount, parentIndex, propInParent;
	NumericMatrix childrenIndexMat,componMappingMat;
  List sepLocation,formula;
	int nDAG;

	depth = as<NumericVector>(treedata["depth"]);
	leaf = as<NumericVector>(treedata["leaf"]);
	nonleaf = as<NumericVector>(treedata["nonleaf"]);
	childrenCount = as<NumericVector>(treedata["childrenCount"]);
	parentIndex = as<NumericVector>(treedata["parentIndex"]);
	propInParent = as<NumericVector>(treedata["propInParent"]);
	childrenIndexMat = as<NumericMatrix>(treedata["childrenIndexMat"]);
  componMappingMat = as<NumericMatrix>(treedata["componMappingMat"]);
	sepLocation = as<List>(treedata["sepLocation"]);
	formula = as<List>(treedata["formula"]);
	nDAG = as<int>(treedata["nDAG"]);


    //if pval is exactly 0 or 1, assign it to a value close to 0 or 1
	int nTree = pval.size();
	for (int i = 0; i < nTree; i++)
	{
		if (pval[i] == 0){
			pval[i] = (1e-11);
		}
		if (pval[i] == 1){
			pval[i] = 1 - (1e-11);
		}
	}

	NumericVector init = init_origin;

	for (int i = 0; i < schedule.size(); i++)
	{

		betaTemp = schedule[i];
		int itertnum;
		if (i == (schedule.size() - 1))
		{
			itertnum = iternum_final;
		}
		else
		{
			itertnum = iternum_hot;
		}
		result = HMTEMC(optim_R, likelihoodfunc_Alt, likelihoodfunc_Null, init, pval, depth, leaf, nonleaf, childrenIndexMat,
			propInParent, childrenCount, parentIndex, nDAG, betaTemp, null_dist, penalty, num_method,
			itertnum, iternum_inner, conv_threshold);
		if ((as<int>(result["iter"]) == iternum_final) && (betaTemp == 1))
		{
			//cout << "EM did not converge in the end!" << endl;
		}
		init = as<NumericVector>(result["par"]);

		//cat("betaTemp:", betaTemp, "; iter:", result$iter, "par:", init, "\n", file=logfile, append=TRUE)

		if (betaTemp < 0.9)
		{
			NumericVector u = runif(9, 0.001, 0.05), temp1(7),temp2(7);
			temp1[0] = u[0];
			temp1[1] = u[1] + 1;
			temp1[2] = u[2];
			temp1[3] = u[3];
			temp1[4] = u[4] + 1;
			temp1[5] = u[5] + 1;
			temp1[6] = 0;
			temp2[0] = 1;
			temp2[1] = 30000;
			temp2[2] = 1 - u[6] * 4;
			temp2[3] = 1;
			temp2[4] = 20 - u[7] * 100;
			temp2[5] = 20 - u[8] * 100;
			temp2[6] = 1 - u[7] * 2;
			init = pmax(init, temp1);
			init = pmin(init, temp2);

		}
	}
	NumericVector cp = as<NumericVector>(result["cp"]);
	NumericVector PDE = getPDEC(cp, parentIndex, componMappingMat, sepLocation, formula);
	return (List::create(Named("par") = result["par"], Named("LL") = result["LL"], Named("iter") = result["iter"], Named("cp") = result["cp"], Named("PDE") = PDE));

}


