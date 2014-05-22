#include<iostream>
#include<primesieve.hpp>
#include<vector>
#include<utility>
#include<cmath>
#include<algorithm>
#include<unordered_map>
#include<set>
#include<list>
using namespace std;
inline long multiplyModuloN(long x,long y,long moduloN)
{
	return ( ( x%moduloN ) * ( y%moduloN ) )%moduloN;
}
inline long addModuloN(long x,long y,long moduloN)
{
	return (  x  +  y) % moduloN;
}
vector<long> generateRanges(long maxNumber,long numRanges)//all ranges have fixed lower end so only upper end is generated
{
long currSlice=1;
vector<long > result(numRanges);
for(unsigned long i=1;i<=numRanges;++i)
		result[i-1]=maxNumber/i;

return result;
}

vector<unordered_map<long,long> >primeFactorize(long maxNumber,const set<long>&primes)
{
	vector<unordered_map<long,long> > result(maxNumber+1);
	unordered_map<long,long> temp;
	bool factorFound=false;
	for( long i=2;i<=maxNumber;++i)
	{
		if(primes.find(i)!=primes.end())
		{
			temp.insert({i,1});
			result[i]=temp;
			temp.clear();
			factorFound=false;
			continue;
		}
		for(auto j=primes.begin();j!=primes.end();++j)
		{
			if( i % (*j)   ==0)
			{
				result[i]=result[i/ (*j)];
				result[i][*j]++;
				factorFound=true;
				break;
			}
			if((*j) > sqrt(i))
				break;
		}
		if(!factorFound)//could not find a factor in primes so number itself is prime
		{
		temp.insert({i,1});
		result[i]=temp;
		temp.clear();
		}
		factorFound=false;
	}
	return result;
}
long intExp(long a,long b,long acc)
{
	if(b==0)
		return acc;
	if(b%2==0)
		return intExp(a*a,b/2,acc);
	else
		return intExp(a,b-1,acc*a);
}
long numRetractionsFunction(long n,const vector<unordered_map<long,long> > & primeFactors)
{
	long result=1;
	for(auto i=primeFactors.at(n).cbegin();i!=primeFactors.at(n).cend();++i)
	{
		result*=(1+intExp(i->first,i->second,1));
	}
	return result-n;
}
inline long signum(long x)
{
	if(x>0)
		return 1;
	if(x<0)
		return -1;
	return 0;
}
void generatePrimeMultipleDivisorsHelp(const unordered_map<long,long> & primeFactors,vector<long>& result,long start,unordered_map<long,long>::const_iterator factorStart)
{
	if (start<result.size())
	{
		for(unsigned long i=0;i<start;++i)
		{
			result[start+i]=-1*(factorStart->first)*result[i];
		}
		generatePrimeMultipleDivisorsHelp(primeFactors,result,2*start,std::next(factorStart));
	}
	 return;
}

vector<long> generatePrimeMultipleDivisors(const unordered_map<long,long>& primeFactors)
{
	vector<long>result(exp2(primeFactors.size()));
	result[0]=1;
	generatePrimeMultipleDivisorsHelp(primeFactors,result,1,primeFactors.begin());
	return result;
}


long calcContributionMultiplier(long maxNumber,long n, const vector<unordered_map<long,long> > & primeFactors,long moduloN)
{
	vector<long> divisors=generatePrimeMultipleDivisors(primeFactors[n]);
	long contrib=0;
	for(auto i=divisors.begin();i<divisors.end();++i)
	{
		contrib+= signum(*i) * ( (maxNumber/(abs(*i)*n) )% moduloN ); ///carefull! Possible source of error if (*i)*n overflows
	}
	return multiplyModuloN(n,(contrib-1),moduloN);
}
long calcSmallNumberContribution(long maxNumber,long maxSmallNumber,const vector<unordered_map<long,long> >& primeFactors,long moduloN)
{
	long contrib=0;
	for(unsigned long i=2;i<=maxSmallNumber;++i)
		contrib=addModuloN(contrib,calcContributionMultiplier(maxNumber,i,primeFactors,moduloN),moduloN);
	return contrib;
}






long rangeContribution(long start,long last,const vector<unordered_map<long,long> > & primeFactors,long maxMultiplier,long moduloN)
{
	if(maxMultiplier==1)
		return 0;
	long total=0;
	// if((start+last)%2==0)
	// 	total=multiplyModuloN(maxMultiplier, multiplyModuloN((start+last)/2,1+last-start,moduloN),moduloN);
	// else
	// 	total=multiplyModuloN(maxMultiplier, multiplyModuloN(start+last,(1+last-start)/2,moduloN),moduloN);
	long startOffset,currStart,lastOffset,currLast,multiplier;
	vector<long> relevantDivisors=generatePrimeMultipleDivisors(primeFactors[maxMultiplier]);
	for(auto iter=relevantDivisors.begin();iter!=relevantDivisors.end();++iter)
	{
		multiplier=abs(*iter);
		if( ( start%multiplier ) != 0 )
		{
			startOffset=(start+multiplier)%multiplier;
			currStart=start+multiplier-startOffset;
		}
		else
			currStart=start;
		
		lastOffset=last%(multiplier);
		currLast=last-lastOffset;
		if(currLast<start)
			continue;//no multiples of this number in the range
		if((currLast+currStart) %2 == 0)
			total += signum(*iter)*	 multiplyModuloN(( currLast+currStart ) / 2,
								   					( currLast - currStart ) / (multiplier )+1,
													 moduloN);

		else
			total += signum(*iter)*multiplyModuloN( currLast+currStart ,
													(( currLast - currStart ) / (multiplier)+1)/2,
													moduloN);
								
		
		
	}
	return total;
	
	
}

long sumRangeContribution(const vector<long>& lastVec,const vector<unordered_map<long,long> >& primeFactors,long moduloN,long rangeMin)
{
	long contrib=0;
	for(unsigned long i=0;i<lastVec.size();++i)
	{
		contrib+= rangeContribution(rangeMin,lastVec[i],primeFactors,i+1,moduloN);
		contrib=contrib%moduloN;
	}
	long leftOverMultiplierMax=lastVec.front()/rangeMin;
	for(unsigned long i=lastVec.size();i<leftOverMultiplierMax;++i)
	{
		contrib+= rangeContribution(rangeMin,rangeMin,primeFactors,i+1,moduloN);
		contrib=contrib%moduloN;	
	}
		
	return contrib;
}

long bruteForceSolution(long maxNumber,long moduloN)
{
	vector<long> primes;
	primesieve::generate_primes(maxNumber,&primes);
	set<long> primeSet(primes.begin(),primes.end());
	auto primeFactors=primeFactorize(maxNumber,primeSet);
	long numRetractions=0;
	for(long i=2;i<=maxNumber;++i)
	{
		numRetractions=addModuloN(numRetractions,numRetractionsFunction(i,primeFactors),moduloN);
	}
	return numRetractions;
}


int main()
{
vector<long> primes;
long moduloN=1000000007;
//long moduloN=83;
long maxNumber;
cin>>maxNumber;
vector<long> ranges=generateRanges(maxNumber,floor(sqrt(maxNumber))-1);

long rangeMinBound=ranges.back()-1;

primesieve::generate_primes(rangeMinBound+1,&primes);

//for(vector<long> ::reverse_iterator i=ranges.rbegin();i!=ranges.rbegin()+10;++i)
//cout<<*i<<" "<< endl;
set<long> primeSet(primes.begin(),primes.end());
auto primeFactors=primeFactorize(rangeMinBound,primeSet);
long largeContrib=sumRangeContribution(ranges,primeFactors,moduloN,rangeMinBound);
long smallContrib=calcSmallNumberContribution(maxNumber,rangeMinBound-1,primeFactors,moduloN);
long ans=addModuloN(addModuloN(largeContrib,smallContrib,moduloN),maxNumber-1,moduloN);// takes into accountcontribution of 1
long test=bruteForceSolution(maxNumber,moduloN);
cout<<ans<<endl;
cout<<test<<endl;
return 0;
}
