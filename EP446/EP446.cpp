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

// implements sieve given here http://www.ams.org/journals/mcom/1959-13-066/S0025-5718-1959-0105784-2/S0025-5718-1959-0105784-2.pdf
inline long multiplyModuloN(long x,long y,long moduloN)
{
	return ( ( x%moduloN ) * ( y%moduloN ) )%moduloN;
}
inline long addModuloN(long x,long y,long moduloN)
{
	return (  x  +  y) % moduloN;
}

unordered_map<long,long> operator+(unordered_map<long,long> x,unordered_map<long,long>y)
{
	unordered_map<long,long> result;
	for(auto i=x.begin();i!=x.end();++i)
			result[i->first]=x[i->first];
	for(auto i=y.begin();i!=y.end();++i)
		result[i->first]+=y[i->first];
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
long intExpModuloN(long a,long b,long acc,long moduloN)
{
	if(b==0)
		return acc;
	if(b%2==0)
		return intExp(multiplyModuloN(a,a,moduloN),b/2,acc);
	else
		return intExp(a,b-1,multiplyModuloN(acc,a,moduloN) );
}
pair<long,long> calcNextABk(long A,long B,long k,long h,long p)
{
	long C=( ( h * ( 1 + A*A ) )/intExp(p,k,1)) %p;
	A=A+C*intExp(p,k,1);
	B=intExp(p,k+1,1)-A;
	return make_pair(A,B);
}
void updateFactorizationMap(unordered_map<long,long>& factors,long prime)
{
	if(factors.find(prime)!=factors.end())
		factors[prime]++;
	else
  	  	factors[prime]=1;
}
vector<unordered_map<long,long> > nSquaredPlusOneSieve(long maxNumber)
{
	vector<long> sieve(maxNumber+1);
	for(unsigned long i=1;i<maxNumber+1;++i)
		sieve[i]=i*i+1;

	vector<unordered_map<long,long> > result(maxNumber+1);
	//divide by 2 where possible
	for(unsigned long i=1;i<maxNumber+1;i+=2)
		{
			if(i<maxNumber+1)
			{
			sieve[i]=sieve[i]/2;
			updateFactorizationMap(result[i],2);
			}
			
		}
	long leftToDo=2;
	long A,B,h,k;
	long currPrime;
	while(leftToDo<maxNumber+1)
	{
		if(sieve[leftToDo]==1)
		{
			leftToDo++;
			continue;
		}
		currPrime=sieve[leftToDo];
			
		A=leftToDo;
		B=currPrime-leftToDo;
		h=( ((currPrime+1)/2)*A )%currPrime;
		k=1;
		while( (A<maxNumber+1) || (B<maxNumber+1) )
		{
			for(unsigned long i=A;i<maxNumber+1;i+=intExp(currPrime,k,1))
			{
				if(i<maxNumber+1)
				{
					updateFactorizationMap(result[i],currPrime);
					if(sieve[i]%currPrime!=0)
						cout<<"error prone sieve-fix it A, i="<<i<<", prime="<<currPrime<< ", k="<<k <<endl;
					sieve[i]/=currPrime;
				}
			}
			for(unsigned long i=B;i<maxNumber+1;i+=intExp(currPrime,k,1))
			{
				if(i<maxNumber+1)
				{
					updateFactorizationMap(result[i],currPrime);
					if(sieve[i]%currPrime!=0)
						cout<<"error prone sieve-fix it B "<<i<<endl;
					sieve[i]/=currPrime;
				}
			}
			if(currPrime > maxNumber+1)
				break;
		auto temp=calcNextABk(A,B,k,h,currPrime);
		A=temp.first;
		B=temp.second;
		k=k+1;

		}
		leftToDo++;
	}
	return result;
}
long numRetractionsFunction(long n,const vector<unordered_map<long,long> > & primeFactors,long moduloN)
{
	long result=1;
	long num=1;
	for(auto i=primeFactors.at(n).cbegin();i!=primeFactors.at(n).cend();++i)
	{
		result=multiplyModuloN(result,(1+intExpModuloN(i->first,i->second,1,moduloN)),moduloN);
		num=multiplyModuloN(num,intExpModuloN(i->first,i->second,1,moduloN),moduloN);
	}
	return result-num;
}
long calcRetractions(long maxNumber,long moduloN)
{
	vector<unordered_map<long,long> > nSquaredPlusOneFactors=nSquaredPlusOneSieve(maxNumber+1);
	vector<unordered_map<long,long> > n4PlusOneFactors(maxNumber+1);
	long result=0;
	n4PlusOneFactors[1]=nSquaredPlusOneFactors[2];
	for(unsigned long  i=1;i<maxNumber;++i)
	{
		n4PlusOneFactors[i+1]=nSquaredPlusOneFactors[i]+nSquaredPlusOneFactors[i+2];
	}

	for(long i=1;i<maxNumber+1;++i)
	{
		result=addModuloN(numRetractionsFunction(i,n4PlusOneFactors,moduloN),result,moduloN);
	}
	return result;
}

int main()
{
	long maxNumber=1e7;
	long moduloN=1000000007;
//	long moduloN=1e4;
	long result=	calcRetractions(maxNumber,moduloN);
	cout<<result<<endl;
}
