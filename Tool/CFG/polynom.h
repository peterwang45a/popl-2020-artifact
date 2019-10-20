#ifndef ____NEUPOLY____
#define ____NEUPOLY____
#define SSP set<set<poly> >
#define polynom poly 
#define term poly_term
#define sort sort_out
#include <iostream>
#include <algorithm>
#include <fstream>
#include <set>
#include <sstream>
#include <vector>
using namespace std;
struct polynom;
struct term //a polynomial term
{
  //friends
  friend polynom substitute(polynom main, string var, polynom rep);
  friend term operator * (term a, term b);
  friend class polynom;
  friend ostream& operator << (ostream &cout, term pterm);
public: //make this private later
  double constant;
  string *variables;
  int *powers;
  int varsize;
  int get_power(string var)
  {
    for(int i=0;i<varsize;++i)
      if(variables[i]==var)
	return powers[i];
    return 0;
  }
  void set_power(string var,int pp)
  {
	bool done=false;
    for(int i=0;i<varsize;++i)
      if(variables[i]==var)
		{powers[i]=pp; done=true;}
	if(!done)
		mult_var(var,pp);
    sort();
  }
  void sort()
  {
    for(int i=0;i<varsize;++i)
      for(int j=i+1;j<varsize;++j)
	if(variables[j]<variables[i] or powers[i]==0)
	  {
	    swap(variables[j],variables[i]);
	    swap(powers[j],powers[i]);
	  }
    //delete all the zeros
    while(varsize!=0 and powers[varsize-1]==0)
      varsize--;
  }
  bool matches(const term &other)
  {
    if(other.varsize!=varsize)
      return false;
    for(int i=0;i<varsize;++i)
      if(other.variables[i]!=variables[i] or other.powers[i]!=powers[i])
	return false;
    return true;
  }
public:
  term(int cst=0) //constructor (with int)
  {
    constant=cst;
    variables=new string[0];
    powers=new int[0];
    varsize=0;
  }
  ~term() //destructor
  {
    delete[] variables;
    delete[] powers;
  }
  void cp_from(const term &t) //copier
  {
    constant=t.constant;
    varsize=t.varsize;
    variables=new string[varsize];
    powers=new int[varsize];
    for(int i=0;i<varsize;++i)
      {
	variables[i]=t.variables[i];
	powers[i]=t.powers[i];
      }
    sort();
  }
  term(const term &t) //copy constructor
  {
    cp_from(t);
  }
  term & operator = (const term &other) //= assignment operator
  {
    delete[] variables;
    delete[] powers;
    cp_from(other);
    return *this;
  }
  bool operator < (const term &other)
  {
    if(varsize!=other.varsize)
      return varsize<other.varsize;
    for(int i=0;i<varsize;++i)
      if(variables[i]!=other.variables[i])
	return variables[i]<other.variables[i];
      else if(powers[i]!=other.powers[i])
	return powers[i]<other.powers[i];
    return false;
  }
  void mult_var(string var,int pw)
  {
    //add var^pw
    for(int i=0;i<varsize;++i)
      if(variables[i]==var)
	{
	  powers[i]+=pw;
	  return;
	}
    string *nvars=new string[varsize+1];
    int *npw=new int[varsize+1];
    for(int i=0;i<varsize;++i)
      {nvars[i]=variables[i]; npw[i]=powers[i];}
    nvars[varsize]=var;
    npw[varsize]=pw;
    varsize++;
    delete[] variables;
    delete[] powers;
    variables=nvars;
    powers=npw;
    sort();
  }
  bool is_constant()
  {
	  return (varsize==0);
  }
};

term operator * (term a, term b)
{
  term ans(a);
  ans.constant*=b.constant;
  for(int i=0;i<b.varsize;++i)
    ans.mult_var(b.variables[i],b.powers[i]);
  ans.sort();
  return ans;
}


struct polynom //a polynomial
{
  friend polynom substitute(polynom main, string var, polynom rep);
  friend polynom operator + (polynom a, polynom b);
  friend polynom operator * (polynom a, term b);
  friend polynom operator * (polynom a, polynom b);
 public: //make private later
  int termsize;
  term *terms;
  void add_term(term t)
  {
    for(int i=0;i<termsize;++i)
      if(terms[i].matches(t))
	{
	  terms[i].constant+=t.constant;
	  return;
	}
    term *novoterms=new term[termsize+1];
    for(int i=0;i<termsize;++i)
      novoterms[i]=terms[i];
    novoterms[termsize]=t;
    termsize++;
    delete[] terms;
    terms=novoterms;
    sort();
  }
public:
 void sort()
  {
    for(int i=0;i<termsize;++i)
      for(int j=i+1;j<termsize;++j)
	if(terms[j]<terms[i] or terms[i].constant==0)
	  swap(terms[j],terms[i]);
    //delete zeros
    while(termsize!=0 and terms[termsize-1].constant==0)
      termsize--;
  }
 int size()
 {
   return termsize;
 }
  ~polynom()//destructor
  {
    delete[] terms;
  }
  polynom()//constructor
  {
    termsize=0;
    terms=new term[0];
  }
  polynom(term t) //constructor with term
  {
    termsize=1;
    terms=new term[1];
    terms[0]=t;
  }
  void cp_from(const polynom &p) //copies p to current item
  {
    termsize=p.termsize;
    terms=new term[termsize];
    for(int i=0;i<termsize;++i)
      terms[i]=p.terms[i];
    sort();
  }
  polynom(const polynom &p) //copy constructor
  {
    cp_from(p);
  }
  polynom & operator = (const polynom &other) //copy assignment operator
  {
    delete[] terms;
    cp_from(other);
    return *this;
  }
  double get_constant()
  {
	  for(int i=0;i<size();++i)
		if(terms[i].is_constant())
			return terms[i].constant;
	  return 0;
  } 
};

polynom operator + (polynom a, polynom b)
{
  polynom ans=a;
  for(int i=0;i<b.termsize;++i)
    ans.add_term(b.terms[i]);
  ans.sort();
  return ans;
}

polynom operator * (polynom a, term b)
{
  polynom ans=a;
  for(int i=0;i<ans.termsize;++i)
    ans.terms[i]=ans.terms[i]*b;
  ans.sort();
  return ans;
}

polynom operator * (polynom a, polynom b)
{
  polynom ans;
  for(int i=0;i<b.termsize;++i)
    ans=ans+(a*b.terms[i]);
  ans.sort();
  return ans;
}

polynom operator - (polynom a)
{
  polynom ans(a);
  ans= ans * (-1);
  ans.sort();
  return ans;
}

polynom operator - (polynom a, polynom b)
{
  polynom ans(a);
  ans=ans+(-b);
  ans.sort();
  return ans;
}

polynom polynom_substitute(polynom main, string var, polynom rep)
{
  polynom ans;
  
  for(int i=0;i<main.termsize;++i)
    {
      polynom add;
      
      term c=main.terms[i];
      int k=0;
      for(int j=0;j<c.varsize;++j)
	if(c.variables[j]==var)
	  {
	    k=c.powers[j];
	    c.powers[j]=0;
	    c.sort();
	  }
      
      add=c;
      while(k--)
	add= add*rep;

      ans=ans+add;
    }
  
  ans.sort();
  return ans;
}


bool operator < (const poly &a, const poly &b)
{
  if(a.termsize!=b.termsize)
    return a.termsize<b.termsize;
  for(int i=0;i<a.termsize;++i)
    if(a.terms[i]<b.terms[i])
      return true;
    else if(b.terms[i]<a.terms[i])
      return false;
  return false;
}

//print operators
ostream& operator << (ostream &cout, term pterm)
{
	//cout<<endl<<" a term:"<<endl;
  cout<<pterm.constant;
  for(int i=0;i<pterm.varsize;++i)
  {
	  if(pterm.variables[i]==" " or pterm.variables[i]=="")
	  {
		cerr<<"ERROR"<<endl;
		exit(0);
		}
    cout<<" * "<<pterm.variables[i]<<"^"<<pterm.powers[i];
  }
  return cout;
}

ostream& operator << (ostream &cout, poly p)
{
  if(p.size()==0)
    cout<<0;
  for(int i=0;i<p.size();++i)
    {
      cout<<p.terms[i];
      if(i+1!=p.size())
	cout<<" + ";
    }
  return cout;
}

poly operator * (double a, poly b)
{
	poly ans;
	term c(a);
	ans=b*c;
	return ans;
}


SSP AND(SSP A, SSP B) //logical and
{
	SSP ANS;
	for(SSP::iterator i=A.begin();i!=A.end();++i)
		for(SSP::iterator j=B.begin();j!=B.end();++j)
		{
			set<poly> a=*i,b=*j;
			a.insert(b.begin(),b.end());
			ANS.insert(a);
		}
	return ANS;
}

bool next_ordering_of(vector<int> &order, const vector<vector<poly> > &vA)
{
	order[order.size()-1]++;
	for(int i=order.size()-1;i>0;--i)
		if(order[i]>=vA[i].size())
		{
			order[i]=0;
			order[i-1]++;
		}
	return (order[0]<vA[0].size());
}

SSP NEG(SSP A) //negation
{
	vector<vector <poly> > vA; //copy it to a vector structrue
	for(SSP::iterator i=A.begin();i!=A.end();++i)
	{
		set<poly> cp=*i;
		vector<poly> cest;
		for(set<poly>::iterator j=cp.begin();j!=cp.end();++j)
			cest.push_back(*j);
		vA.push_back(cest);
	}
	
	SSP ANS;
	//create the negation
	vector<int> order(vA.size(),0);
	do
	{
		set<poly> sp;
		for(int i=0;i<vA.size();++i)
			sp.insert((-1)*vA[i][order[i]]);
		ANS.insert(sp);
	}while(next_ordering_of(order,vA));
	return ANS;
}


poly read_poly(string source) //read one polynomial from the string
{
  source=source+" + "; //add a plus at the end
  stringstream cin;
  cin<<source;
  poly ans;
  poly_term cterm;
  string token,prevtoken="+";
  while(cin>>token)
    {
      if(token=="+")
	ans=ans+cterm;
      else if(prevtoken=="+")
	{
	  //this is a constant
	  stringstream cnv;
	  cnv<<token;
	  poly_term novo; cterm=novo; //clear cterm
	  cnv>>cterm.constant; //read the constant
	}
      else if(prevtoken=="*")
	{
	  //it's a variable and its power
	  for(int i=0;i<token.length();++i)
	    if(token[i]=='^')
	      token[i]=' ';
	  stringstream cnv;
	  cnv<<token;
	  string variable; int power;
	  cnv>>variable>>power;
	  cterm.mult_var(variable,power);
	}

      prevtoken=token;
    }
  return ans;
}


#undef polynom
#undef sort
#undef term
#endif
