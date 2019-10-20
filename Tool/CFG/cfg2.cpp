//Parser
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <vector>
#include <string>
#include <sstream>
#include <cctype>
#include <set>
#include <map>
#include <fstream>
#include <deque>
#include <stack>
#include <unordered_map>
#include <tuple>
#include <list>
#include <cstring>
#include "polynom.h"
using namespace std;

#define vector deque

#define MAXL 100000 //Maximum length of the program
char input[MAXL];
string program;

#define part(x,a,b) (x.substr((a),((b)-(a))))
#define pb push_back
#define make_edge(a,b,c) edge((int *)(a),(int *)(b),(c))
#define mp make_pair
const int MAXD=2; //maximum degree

int branchnum=1;
int assnum=0;


set<string> pvars; //program variables
map<int, pair<double, double> > rvars; //randomized variables
map<void *, int> label_map; //labels
//New add
//map<string, int> gamma_number;
map<string,string> gamma_matrix;
int max_label = 0;

//ofstream fout("parser_output.txt"); //output
ofstream fout;
void print_pair(const set<set<poly> > &set_set, poly f) //prints pairs to output
{
	fout<<"NEW PAIR:"<<endl;
	for(set<set<poly> >::iterator iter1=set_set.begin();iter1!=set_set.end();++iter1)
	{
		for(set<poly>::iterator iter=iter1->begin();iter!=iter1->end();++iter)
			fout<<*iter<<endl;
	}
	fout<<"@@@@@@@@@@"<<endl;
	fout<<f<<endl;
	fout<<"@@@@@@@@@@"<<endl;
}

string ascize(string x)
{
	stringstream cnv;
	for(int i=0;i<x.length();++i)
		cnv<<int(x[i])<<"--";
	cnv>>x;
	return x;
}

void substitute_random_variables_with_values(poly_term &term) //works for at-most-quadratic polynomials only
{
	for(map<int, pair<double, double> >::iterator iter=rvars.begin();iter!=rvars.end();++iter)
	{
		string var; //variable name
		stringstream cnv;
		cnv<<"r_"<<iter->first;
		cnv>>var;
		
		double a=iter->second.first,b=iter->second.second;
		if(term.get_power(var)==1)
		{
		  term.set_power(var,0);
			term.constant*=(a+b)/2;
		}
		if(term.get_power(var)==2)
		{
		  term.set_power(var,0);
			double ex=(a+b)/2;
			double variance=(b-a)*(b-a)/12;
			term.constant*=variance + ex * ex ; //E[X^2] = Var(X) + E[X]^2
		}
	}
	
	term.sort_out();
}

void substitute_random_variables_with_values(poly &p) //works for at-most-quadratic polynomials only
{
	for(int termi=0;termi<p.size();++termi)
	{
		poly_term &term=p.terms[termi];
		substitute_random_variables_with_values(term);
	}
	p.sort_out();
}

double to_double(string sth)
{
	stringstream cnv;
	cnv<<sth;
	double ret;
	cnv>>ret;
	return ret;
}

void add_pvar(string p)
{
	/*while(!isdigit(p[0]))
		p=p.substr(1,p.length());
	stringstream cnv;
	cnv<<p;
	int x;
	cnv>>x;*/
	stringstream cnv;
	cnv<<p;
	cnv>>p; //This gets rid of the spaces
	pvars.insert(p);
	//gamma_number[p]=0;
	gamma_matrix[p]="";
}

void add_rvar(string r)
{
	for(int i=0;i<r.length();++i)
		if(!isdigit(r[i]) and r[i]!='.')
			r[i]=' ';
	int a; double b=0,c=0;
	stringstream cnv;
	cnv<<r;
	cnv>>a;
	if(cnv>>b>>c)
	{
		rvars[a]=mp(b,c);
	}
	else
	{
		b=c=0;
		if(rvars.find(a)==rvars.end())
			rvars[a]=mp(b,c);
    }
}

struct edge //CFG edge
{
  edge(int *x,int *y, string t)
  {
    u=x;
    v=y;
    text=t;
  }
  int*u,*v;
  string text;
};

vector<edge> CFG;

//New Add
//Struct for storing infomation
struct info
{
  string type; //type of distribution
  string parameters;
  string evaluations;
  //constructor
  info()
  {
	type="";
	parameters="";
	evaluations="";
  }
  info(string t,string p,string e)
  {
	type=t;
	parameters=p;
	evaluations=e;
  }
  void setType(string t)
  {
  	type=t;
  }
  void setPara(string p)
  {
  	parameters=p;
  }
  void setEval(string e)
  {
  	evaluations=e;
  }
  void addPara(string p)
  {
  	parameters+=" "+p;
  }
  void addEval(string e)
  {
  	evaluations+=" "+e;
  }
  string getType()
  {
  	return type;
  }
  string getPara()
  {
  	return parameters;
  }
  string getEval()
  {
  	return evaluations;
  }
  string printInfo()
  {
  	string str="";
  	if(evaluations!="")
  		str=str+evaluations+"\n";
  	str=str+parameters+"\n";
  	return str;
  }
};

unordered_map<string,info> distributions;

void skip_spaces(int &begin, int &end) //skips the spaces in the beginning and the end
{
  while(begin<program.size() and isspace(program[begin]))
    begin++;
  while(end>0 and isspace(program[end-1]))
    end--;
}

bool is_good_paren(int begin, int end) //whether this has a valid parentheses combination
{
  int cnt=0;
  for(int i=begin;i<end;++i)
    if(program[i]=='(')
      cnt++;
    else if(program[i]==')')
      {
	cnt--;
	if(cnt<0)
	  return false;
      }
  if(cnt==0)
    return true;
  return false;
}

bool next_combination_with_this_sum(int *begin,int size)
{
	if(size==1)
		return false;
	if(next_combination_with_this_sum(begin+1,size-1))
		return true;
	if((*begin)==0)
		return false;
	(*begin)--;
	(*(begin+1))++;
	for(int i=2;i<size;++i)
	{
		(*(begin+1))+=(*(begin+i));
		(*(begin+i))=0;
	}
	return true;
}


int last_used_label=0;
struct node
{
  int label; //label
  string type; //The type of entity that this node corresponds to
  string constant; //The constant factor (subtype) of this node
  int begin,end; //begin and end point in program
  vector<node *> children; //children in the parse tree
  node *bracket; //bexpr written in bracket next to this node (for stmt only)
  node *next_statement; //for statements only
  poly polynom; //corresponding polynomial -- for expr and rexpr
  set<poly> polynom_set; //correspoding set of polynomils (positivity) -- for polyexpr
  set<set<poly> > polynom_set_set; //corresponding set of sets of polynomials (positivity) -- for bexpr
  poly monomials; //contains all the monomials corresponding to the current label
  #define eta monomials
  poly pre_eta; //pre_eta as in the paper
  //constructor
  node(string t)
  {
	  bracket=NULL;
	  begin=end=-1;
	  type=t;
	  constant="";
	  label=0;
  }
  node(string t, int b, int e, node *nx=NULL)
  {
	label=0;
    bracket=NULL;
    next_statement=nx;
    type=t;
    begin=b;
    end=e;
    skip_spaces(begin,end);
    if(begin>end)
      cerr<<"Error! begin= "<<b<<" end= "<<e<<endl;
    process();
  }
  
  //recursive functions
  void recursively_create_monomials()
	{
		if(pvars.size()==0)
			return;
		if(label!=0)
		{
			int total_counter=0;
			//create the monomials here
			for(int cnt=0;cnt<=MAXD;++cnt)
			{
				int put[pvars.size()];
				fill(put,put+pvars.size(),0);
				put[0]=cnt;
				do
				{
					poly_term add;
					add.constant=1;
					//add the variables
					set<string>::iterator iter=pvars.begin();
					int I=0;
					for(;iter!=pvars.end();I++,++iter)
					{
						//add.variables[*iter]+=put[I];
						add.set_power(*iter,put[I]+add.get_power(*iter));
					}
					//add the fresh variable (constant factor)
					stringstream cnv;
					cnv<<"a_"<<label<<","<<total_counter;
					string v;
					cnv>>v;
					//add.variables[v]++;
					add.set_power(v,add.get_power(v)+1);
					total_counter++;
					
					//add this term to the polynomial
					add.sort_out();
					//monomials.push_back(add);
					monomials=monomials+add;
					monomials.sort_out();
					
				}while(next_combination_with_this_sum(put,pvars.size()));
				
			}
			
			cout<<"Label: "<<label<<"  monomials:"<<endl;
			cout<<monomials<<endl;
			
		}
		//recurse
		for(int i=0;i<children.size();++i)
			children[i]->recursively_create_monomials();
	}
	
	void recursively_calculate_pre_etas()
	{
		bool calc=false;
		//check if this is a probabilistic if
		if(type=="stmt" and constant=="if" and children[0]->constant=="prob")
		{
			calc=true;
			for(int i=0;i<CFG.size();++i)
				if(((node *)CFG[i].u)==this)
				{
					node *next=(node *)(CFG[i].v);
					double pr;
					stringstream cnv; string e=CFG[i].text;
					for(int j=0;j<e.length();++j)
						if(!isdigit(e[j]) and e[j]!='.' and e[j]!='/')
							e[j]=' ';
					cnv<<e;
					cnv>>pr;
					if(next!=NULL)
					pre_eta=pre_eta+pr*(next->eta);
				}
		}
		if(type=="stmt" and constant==":=")
		{
			//set pre_eta to eta of the next statement according to CFG
			for(int i=0;i<CFG.size();++i)
				if(((node *)CFG[i].u)==this)
				{
					node *next=(node *)(CFG[i].v);
					if(next!=NULL)
					pre_eta=next->eta;
				}
			
			//calculate left and right hand sides
			calc=true;
			vector<string> left_hand_side;
			vector<poly> right_hand_side;
			//create the left hand side
			node *pointer=children[0];
			while(pointer->type=="pvarlist")
			{
				left_hand_side.pb(pointer->children[0]->constant);
				if(pointer->children.size()>1)
				pointer=pointer->children[1];
				else break;
			}
			
			//create the right hand side
			pointer=children[1];
			while(pointer->type=="rexprlist")
			{
				right_hand_side.pb(pointer->children[0]->polynom);
				if(pointer->children.size()>1)
				pointer=pointer->children[1];
				else break;
			} 
		
		
			
			/*//print it
			cout<<"left hand side:"<<endl;
			for(int i=0;i<left_hand_side.size();++i)
				cout<<left_hand_side[i]<<endl;
				
			cout<<"right hand side:"<<endl;
			for(int i=0;i<right_hand_side.size();++i)
				cout<<right_hand_side[i]<<endl;*/
				
			//subsitute random variables in the right hand side with their values
			for(int i=0;i<right_hand_side.size();++i)
				substitute_random_variables_with_values(right_hand_side[i]);
				
			//substitute left hand side by right hand side in the pre_eta
			for(int i=0;i<left_hand_side.size() and i<right_hand_side.size();++i)
				pre_eta=polynom_substitute(pre_eta, left_hand_side[i], right_hand_side[i]);
				
			
		}
		if(calc) //print it
			cout<<"pre_eta ("<<label<<"): "<<pre_eta<<endl;
		//recurse
		for(int i=0;i<children.size();++i)
			children[i]->recursively_calculate_pre_etas();
	}
	
private:
  //procs
  void proc_pvar()
  {
    while(begin<program.length() and program[begin]!='p')
      begin++;
  // Modified by QXD
    while(end>0 and !isalnum(program[end-1]))
      end--;
    constant=part(program,begin,end);
    add_pvar(part(program,begin,end));
  }

  void proc_rvar()
  {
    while(begin<program.length() and program[begin]!='r')
      begin++;
  // Modified by QXD
    while(end>0 and !isalnum(program[end-1]))
      end--;
    constant=part(program,begin,end);
    add_rvar(part(program,begin,end));
  }

  void proc_constant()
  {
    while(begin<program.length() and !isdigit(program[begin]) and program[begin]!='.')
      begin++;
    while(end>0 and !isdigit(program[end-1]))
      end--;
    constant=part(program,begin,end);
    if(constant[0]=='.')
      constant="0"+constant;
  }

  void proc_pvarlist()
  {
    int comma=-1;
    for(int i=begin;i<end;++i)
      if(program[i]==',')
	{comma=i; break;}
    if(comma==-1)
      {
	constant="single pvar";
	node *ch = new node("pvar",begin,end);
	children.pb(ch);
      }
    else
      {
	constant = "multiple pvar";
	node *ch= new node("pvar",begin,comma);
	children.pb(ch);
	node *ch2= new node("pvarlist",comma+1,end);
	children.pb(ch2);
      }
  }

  void proc_rexprlist()
  {
    int comma=-1;
    for(int i=begin;i<end;++i)
      if(program[i]==',')
	{comma=i; break;}
    if(comma==-1)
      {
	constant = "single rexpr";
	node *ch=new node("rexpr",begin,end);
	children.pb(ch);
      }
    else
      {
	constant = "multiple rexpr";
	node *ch= new node("rexpr", begin, comma);
	children.pb(ch);
	node *ch2=new node("rexprlist",comma+1,end);
	children.pb(ch2);
      }
  }

  void expr_create_polynom(bool rex=false) //for both expr and rexpr
  {
	  if(constant=="()")
	  {
		  polynom=children[0]->polynom;
		  return;
	  }
	  if(constant=="*")
	  {
		  polynom=children[0]->polynom*children[1]->polynom;
		  return;
	  }
	  if(constant=="+")
	  {
		  polynom=children[0]->polynom+children[1]->polynom;
		  return;
	  }
	  if(constant=="-")
	  {
		  polynom=children[0]->polynom+(-1)*children[1]->polynom;
		  return;
	  }





	  if(constant=="single pvar")
	  {
		  string pvarname=children[0]->constant;
		  poly_term term;
		  term.constant=1;
		  //term.variables[pvarname]=1;
		  term.set_power(pvarname,1);
		  polynom=polynom+term;
		  return;
	  }
	  if(constant=="single rvar")
	  {
		  string rvarname=children[0]->constant;
		  for(int i=0;i<rvarname.size();++i)
			if(rvarname[i]=='{')
				rvarname=rvarname.substr(0,i);
		  poly_term term;
		  term.constant=1;
		  //term.variables[rvarname]=1;
		  term.set_power(rvarname,1);
		  //polynom.pb(term);
		  polynom=polynom+term;
		  return;
	  }
	  if(constant=="single constant")
	  {
		  string cst=children[0]->constant;
		  poly_term term;
		  term.constant=to_double(cst);
		  //polynom.pb(term);
		  polynom=polynom+term;
		  return;
	  }
  }

  //Modified by QXD
  void proc_expr(bool allow_rval=0) //for both expr and rexpr
  {
    skip_spaces(begin,end);
    //check for parentheses
    if(begin<program.size() and program[begin]=='(' and end>0 and program[end-1]==')' and is_good_paren(begin+1,end-1))
      {
	constant="()";
	node *ch=new node(type,begin+1,end-1);
	children.pb(ch);
	return;
      }
    //check for +,-
    int cntopen=0; //number of open parentheses
    int plusminus=-1; //place of + or - sign
    int mult=-1; //place of * sign
    for(int i=begin;i<end;++i)
      if(program[i]=='(')
	cntopen++;
      else if(program[i]==')')
	{
	  cntopen--;
	  if(cntopen<0)
	    cerr<<"Error: Cannot find '(' for ')' at position:"<<i<<endl;
	}
      else if((program[i]=='+' or program[i]=='-') and cntopen==0)
	//{
		plusminus=i; 
	//	break;}
      else if(program[i]=='*' and cntopen==0)
	mult=i;
    
    if(plusminus!=-1)
      {
	constant=part(program,plusminus,plusminus+1);
	children.resize(2);
	children[0]=new node(type,begin,plusminus);
	children[1]=new node(type,plusminus+1,end);
	return;
      }

    if(mult!=-1)
      {
	constant="*";
	children.resize(2);
	children[0]=new node(type,begin,mult);
	children[1]=new node(type,mult+1,end);
	return;
      }

    if(program[begin]=='p')
      {
	constant = "single pvar";
	node *ch=new node("pvar",begin,end);
	children.pb(ch);
      }
    else if(program[begin]=='r' and allow_rval)
      {
	constant = "single rvar";
	node *ch= new node("rvar",begin,end);
	children.pb(ch);
      }
    else
      {
	constant = "single constant";
	node *ch=new node("constant",begin,end);
	children.pb(ch);
      }

 }//proc_expr

  void proc_literal()
  {
    int sign=-1;
    for(int i=begin;i<end;++i)
      if(part(program,i,i+2)=="<=" or part(program,i,i+2)==">=")
	{sign = i; break;}
    if(sign!=-1)
      {
	constant=part(program,sign,sign+2);
	children.resize(2);
	children[0]=new node("expr",begin,sign);
	children[1]=new node("expr",sign+2,end);
      }
  }
 

  void normalize_literal() //convert to the format sth>=0
  {
	//create a zero expr
    node *exprzero=new node("expr");
    exprzero->constant="single constant";
    node *ch=new node("constant");
    ch->constant="0";
    exprzero->children.pb(ch);
    exprzero->expr_create_polynom();
    
    //create a minus expr
    node *minus=new node("expr");
    minus->constant="-";
    minus->children.resize(2);
    minus->children[0]=(constant=="<=")?children[1]:children[0];
    minus->children[1]=(constant=="<=")?children[0]:children[1];
    minus->expr_create_polynom();
    
    constant=">=";
    children[0]=minus;
    children[1]=exprzero;
  }

  void proc_polyexpr()
  {
    int AND=-1;
    for(int i=begin;i<end;++i)
      if(part(program,i,i+3)=="and")
	{AND=i; break;}
    if(AND==-1)
      {
	constant = "single literal";
	node *ch=new node("literal",begin,end);
	children.pb(ch);
      }
    else
      {
	constant = "and";
	children.resize(2);
	children[0]=new node("literal",begin,AND);
	children[1]=new node("polyexpr",AND+3,end);
      }
  }
  
  void polyexpr_create_polynomial_set()
  {
	  for(int i=0;i<children.size();++i)
		if(children[i]->type=="literal")
		{
			poly add=children[i]->children[0]->polynom;
			//sort(add.begin(),add.end());
			add.sort_out();
			polynom_set.insert(add);
		}
		else
			for(set<poly>::iterator iter=children[i]->polynom_set.begin();iter!=children[i]->polynom_set.end();++iter)
				polynom_set.insert(*iter);
  }

  void proc_bexpr()
  {
    int OR=-1;
    for(int i=begin;i<end;++i)
      if(part(program,i,i+2)=="or")
	{OR=i; break;}
    if(OR==-1)
      {
	constant = "single polyexpr";
	node *ch=new node("polyexpr",begin,end);
	children.pb(ch);
      }
    else
      {
	constant = "or";
	children.resize(2);
	children[0]=new node("bexpr",begin,OR);
	children[1]=new node("bexpr",OR+2,end);
      }
  }
  
  void bexpr_create_polynomial_set_of_set()
  {
	  for(int i=0;i<children.size();++i)
		if(children[i]->type=="polyexpr")
			polynom_set_set.insert(children[i]->polynom_set);
		else
			for(set<set<poly> >::iterator iter=children[i]->polynom_set_set.begin();iter!=children[i]->polynom_set_set.end();++iter)
				polynom_set_set.insert(*iter);
  }
  

  void proc_ndbexpr()
  {
    skip_spaces(begin,end);
    if(part(program,begin,end)=="demon")
      {
	constant="demon";
	return;
      }
    if(part(program,begin,begin+4)=="prob")
      {
	constant="prob";
	int open=-1,close=-1;
	for(int i=begin;i<end;++i)
	  if(program[i]=='(')
	    open=i;
	  else if(program[i]==')')
	    close=i;
	if(open==-1 or close==-1)
	  cerr<<"Invalid probability"<<endl;
	node *ch=new node("constant",open+1,close);
	children.pb(ch);
	return;
      }
    constant="single bexpr";
    children.resize(1);
    children[0]=new node("bexpr",begin,end);
  }

  //New add
  void supplement_missing_next_branching(node *first_branch,node *second_branch,node *next,int next_label)
  {
	while(first_branch!=NULL and first_branch->constant=="several statements")
		first_branch=first_branch->children[1];
	while(second_branch!=NULL and second_branch->constant=="several statements")
		second_branch=second_branch->children[1];
  	//cout<<first_branch->label<<" "<<second_branch->label<<endl;
  	if(first_branch->constant=="if")
  	{
  		node *f=first_branch->children[1];
		node *s=first_branch->children[2];
  		supplement_missing_next_branching(f,s,next,next_label);
  	}
  	else
  	{
  		int current_label1=first_branch->label;
	  	for(int i=0;i<CFG.size();++i)
		{
			//cout<<"========"<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
			//cout<<"========"<<CFG[i].u<<" -> "<<CFG[i].v<<" :: "<<CFG[i].text<<endl;
			if(label_map[CFG[i].u]==current_label1 && label_map[CFG[i].v]==0)
			{
				CFG[i].v=(int *)next;
				label_map[CFG[i].v]=next_label;
			}
			//cout<<"--------"<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
			//cout<<"--------"<<CFG[i].u<<" -> "<<CFG[i].v<<" :: "<<CFG[i].text<<endl;
		}
  	}
  	if(second_branch->constant=="if")
  	{
  		node *f=second_branch->children[1];
		node *s=second_branch->children[2];
  		supplement_missing_next_branching(f,s,next,next_label);
  	}
  	else
  	{
  		int current_label2=second_branch->label;
	  	for(int i=0;i<CFG.size();++i)
		{
			//cout<<"========"<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
			//cout<<"========"<<CFG[i].u<<" -> "<<CFG[i].v<<" :: "<<CFG[i].text<<endl;
			if(label_map[CFG[i].u]==current_label2 && label_map[CFG[i].v]==0)
			{
				CFG[i].v=(int *)next;
				label_map[CFG[i].v]=next_label;
			}
			//cout<<"--------"<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
			//cout<<"--------"<<CFG[i].u<<" -> "<<CFG[i].v<<" :: "<<CFG[i].text<<endl;
		}
  	}
  }
  
  //New add
  void supplement_missing_next(node *child)
  {
  	  node *next = child->next_statement;
	  int next_label=child->next_statement->label;
	  if(child->constant=="if")
	  {
	  	  node *first_branch=child->children[1];
		  node *second_branch=child->children[2];
	  	  supplement_missing_next_branching(first_branch,second_branch,next,next_label);
	  }
	  else
	  {
	  	  int current_label=child->label;
		  for(int i=0;i<CFG.size();++i)
		  {
		  	if(label_map[CFG[i].u]==current_label && label_map[CFG[i].v]==0)
		  	{
		  		CFG[i].v=(int *)child->next_statement;
		  		label_map[CFG[i].v]=next_label;
		  	}
		  }
	  }
	  return;
  }

  void proc_stmt()
  {
    skip_spaces(begin,end);








    //look for semicolons
    int open=0; //open if's and while's
    for(int i=begin;i<end;++i)
      if(part(program,i,i+2)=="if")
	open++;
      else if(part(program,i,i+5)=="while")
	open++;
      else if(part(program,i,i+2)=="od" or part(program,i,i+2)=="fi")
	open--;
      else if(open==0 and program[i]==';')
	{
	  constant = "several statements";
	  children.resize(2);
	  //New add
	  //Switch label into a correct order
	  children[0]=new node("stmt",begin,i,children[1]);
	  children[1]=new node("stmt",i+1,end,next_statement);
	  children[0]->next_statement=children[1];
	  //New add
	    for(int i=0;i<CFG.size();++i)
	    {
	      //cout<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
	      edge &e=CFG[i];
	      node *u=(node *)e.u,*v=(node *)e.v;
	      //cout<<u<<endl;
	      //cout<<v<<endl;
	      while(u!=NULL and u->constant=="several statements")
			u=u->children[1];
	      while(v!=NULL and v->constant=="several statements")
			v=v->children[0];
	      e.u = (int *)u;
	      e.v = (int *)v;
		  //cout<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
	    }
	  supplement_missing_next(children[0]);
	  return;
	}

    //look for brackets
    skip_spaces(begin,end);
    if(program[begin]=='[')
      {
	int closed_bracket=-1;
	for(int j=begin+1;j<end;++j)
	  if(program[j]==']')
	    {closed_bracket=j; break;}
	bracket=new node("bexpr",begin+1,closed_bracket);
	begin=closed_bracket+1;
      }
    skip_spaces(begin,end);

    //check if it is a while
    if(part(program,begin,begin+5)=="while")
      {
	int firstdo=-1,lastod=-1;
	for(int i=begin;i<end;++i)
	  if(part(program,i,i+2)=="do")
	    {firstdo=i; break;}
	for(int i=end-2;i>=begin;--i)
	  if(part(program,i,i+2)=="od")
	    {lastod=i; break;}
	if(firstdo==-1 or lastod==-1)
	  cerr<<"Bad while in range ["<<begin<<","<<end<<")"<<endl;
	constant = "while";
	label=++last_used_label;
	label_map[this]=label;
	children.resize(2);
	children[0]=new node("bexpr",begin+5,firstdo);
	children[1]=new node("stmt",firstdo+2,lastod,this);
	CFG.pb(make_edge(this,children[1],part(program,children[0]->begin,children[0]->end)));
	CFG.pb(make_edge(this,next_statement,"~("+part(program,children[0]->begin,children[0]->end)+")"));
	return;
      }

    //check if it is an if
    if(part(program,begin,begin+2)=="if")
      {
    int star=0;
	int firstthen=-1; //first then
	for(int i=begin;i<end;++i)
	  if(part(program,i,i+4)=="then")
	    {firstthen=i; 
			for(int ii=begin;ii<i+1;++ii)
			  if(part(program,ii,ii+1)=="*")
			    star=1;
		break;}
	int lastfi=-1;
	for(int i=end-1;i>=begin;--i)
	  if(part(program,i,i+2)=="fi")
	    {lastfi=i; break;}
	if(lastfi==-1 or firstthen==-1)
	  cerr<<"Invalid if between ["<<begin<<", "<<end<<")"<<endl;
	int correselse=-1; //corresponding else
	int ifcnt=0;
	for(int i=begin;i<end;++i)
	  if(part(program,i,i+2)=="if")
	    ifcnt++;
	  else if(part(program,i,i+4)=="else")
	    {
	      ifcnt--;
	      if(ifcnt==0)
		{
		  correselse=i;
		  break;
		}
	    }
	constant="if";
	label=++last_used_label;
	label_map[this]=label;
	children.resize(3);
	if(star == 0)
	{
		children[0]=new node("ndbexpr",begin+2,firstthen);
		children[1]=new node("stmt",firstthen+4,correselse,next_statement);
		children[2]=new node("stmt",correselse+4,lastfi,next_statement);
		CFG.pb(make_edge(this,children[1],part(program,children[0]->begin,children[0]->end)));
		CFG.pb(make_edge(this,children[2],"~("+part(program,children[0]->begin,children[0]->end)+")"));
	}
	else
	{
		children[0]=new node("star");
		children[1]=new node("stmt",firstthen+4,correselse,next_statement);
		children[2]=new node("stmt",correselse+4,lastfi,next_statement);
		CFG.pb(make_edge(this,children[1],"(*)"));
		CFG.pb(make_edge(this,children[2],"~(*)"));
	}
	return;
      }

    //New add
    //tick(rexpr)
    if(part(program,begin,begin+4)=="tick")
      {
    int open=0;
    int close=0;
    for(int i=begin;i<end;++i)
	  if(program[i]=='(')
	    open=i;
	  else if(program[i]==')')
	    close=i;
	if(open==-1 or close==-1)
	  cerr<<"Invalid probability"<<endl;
    constant="tick";
	label=++last_used_label;
	label_map[this]=label;
	children.resize(1);
    children[0]=new node("rexpr",open+1,close);
    CFG.pb(make_edge(this,next_statement,"true"));
	return;
       }

    //New Modify
    //Moved by QXD
    if(part(program,begin,begin+4)=="skip")
      {
	constant="skip";
	label=++last_used_label;
	label_map[this]=label;
	CFG.pb(make_edge(this,next_statement,"true"));
	return;
      }

    //New add
    //end
    if(part(program,begin,begin+3)=="end")
      {
    constant="end";
	label=++last_used_label;
	label_map[this]=label;
	return;
       }

    //pvarlist := rexprlist
    int eq=-1;
    for(int i=begin;i<end;++i)
      if(part(program,i,i+2)==":=")
	{eq=i; break;}
    if(eq==-1)
      cerr<<"Invalid stmt in ["<<begin<<", "<<end<<")"<<endl;
    constant=":=";
    label=++last_used_label;
    label_map[this]=label;
    children.resize(2);
    children[0]=new node("pvarlist",begin,eq);
    children[1]=new node("rexprlist",eq+2,end);
    CFG.pb(make_edge(this,next_statement,"true"));
    return;
  }

  //main process
  void process()
  {
    if(type=="pvar")
      proc_pvar();
    else if(type=="rvar")
      proc_rvar();
    else if(type=="constant")
      proc_constant();
    else if(type=="pvarlist")
      proc_pvarlist();
    else if(type=="expr")
    {
      proc_expr();
      expr_create_polynom();
	}
    else if(type=="rexpr")
    {
      proc_expr(1);
      expr_create_polynom(1);
    }
    else if(type=="rexprlist")
      proc_rexprlist();
    else if(type=="literal")
      {
	proc_literal();
	normalize_literal();
      }
    else if(type=="polyexpr")
    {
      proc_polyexpr();
      polyexpr_create_polynomial_set();
    }
    else if(type=="bexpr")
    {
      proc_bexpr();
      bexpr_create_polynomial_set_of_set();
	}
    else if(type=="ndbexpr")
      proc_ndbexpr();
    else if(type=="stmt")
      proc_stmt();
    else
      cerr<<"Undefined type between "<<begin<<" "<<end<<endl;
  }
public:


void recursively_write_output()
{
	//define epsilon
	poly_term epsilon_term;
	epsilon_term.constant=1;
	poly epsilon;
	//epsilon.pb(epsilon_term);
	epsilon=epsilon+epsilon_term;
	
	//part 2
	if(type=="stmt" and constant=="if" and children[0]->constant=="prob") //probabilistic if
		{
			poly f;
			f=eta;
			f=f+(-1)*pre_eta;
			f=f+epsilon;
			print_pair(bracket->polynom_set_set,f);
		}
	
	//part 3
	if(type=="stmt" and constant==":=")
		print_pair(bracket->polynom_set_set,eta+((-1)*pre_eta)+epsilon);
		
	//part 4
	if(type=="stmt" and constant=="if" and children[0]->constant=="single bexpr")
	{
		//make I and phi
		SSP I = bracket->polynom_set_set;
		SSP phi = children[0]->children[0]->polynom_set_set;
		
		//make (I1, F1)
		SSP I1=AND(I,phi);
		poly F1=eta + ((-1)*(children[1]->eta)) + epsilon;
		
		//make (I2,F2)
		SSP I2=AND(I,NEG(phi));
		poly F2=eta + ((-1)*(children[2]->eta)) + epsilon;
		
		print_pair(I1,F1);
		print_pair(I2,F2);
	}
	
	//part 5
	if(type=="stmt" and constant=="if" and children[0]->constant=="demon")
	{
		SSP I = bracket->polynom_set_set;
		poly F1,F2;
		F1= children[1]->eta + ((-1)*eta) + epsilon;
		F2= children[2]->eta + ((-1)*eta) + epsilon;
		print_pair(I,F1);
		print_pair(I,F2);
	}	
	
	//part 6
	if(label!=0)
		print_pair(bracket->polynom_set_set, eta);
	
	//recurse
	for(int i=0;i<children.size();++i)
		children[i]->recursively_write_output();
}

void set_set_print(ostream &wr = cout)
{
	wr<<"Polynom_set_set:"<<endl;
		for(set<set<poly> >::iterator iter1=polynom_set_set.begin();iter1!=polynom_set_set.end();++iter1)
		{
			for(set<poly>::iterator iter=iter1->begin();iter!=iter1->end();++iter)
				wr<<"|"<<*iter<<"|"<<endl;
		}
}

  void print()
  {
    if(label!=0)
		cout<<"Label: "<<label<<endl;
	if(next_statement!=NULL)
		cout<<"Next Label: "<<next_statement->label<<endl;
    cout<<"Add: "<<this<<"\t";
    cout<<"Type: "<<type<<"\t";
    cout<<"Range: ["<<begin<<", "<<end<<")\t";
    cout<<"Const: |"<<constant<<"|"<<endl;
    // New add
    if(type=="stmt")
    {
    	if(constant=="while")
    	{
    		//fout<<"branching"<<endl;
    		// for(int i=0;i<CFG.size();++i)
    		// {
    		// 	if(label_map[CFG[i].u]==label)
    		// 	{
    		// 		if(label_map[CFG[i].v]!=0)
    		// 		{
    		// 			fout<<label_map[CFG[i].v]<<" ";
    		// 		}
    		// 		else
    		// 		{
    		// 			cout<<"Error: "<<label<<" misses a next label."<<endl;
    		// 		}
    		// 	}
    		// }
      // 		fout<<endl;
      		//New add
			//Gamma
   //  		string r="";
			// set<string>::iterator iter=pvars.begin();
		 //    for(;iter!=pvars.end();++iter)
		 //    {
		 //    	gamma_matrix[*iter].clear();
		 //    }
		 //    //cout<<"Run extract_gamma_coeffs"<<endl;
		 //    string res="";
			// int sum = 0;
		 //    res+=extract_gamma_coeffs(bracket,1,"",sum);
		 //    r=extract_gamma_coeffs(children[0],1,res,sum);
		 //    if(r!="")
		 //    	res+=r;
			// fout<<sum<<endl;
		 //    fout<<res;
			// //fout<<"[#Input current value of pv there.]"<<endl;
			// iter=pvars.begin();
		 //    for(;iter!=pvars.end();++iter)
		 //    {
		 //    	gamma_matrix[*iter].clear();
		 //    }
		 //    //cout<<"Run extract_gamma_coeffs"<<endl;
		 //    res.clear();
			// sum=0;
			// res+=extract_gamma_coeffs(bracket,1,"",sum);
		 //    node *medi=next_statement;
		 //    while(medi!=NULL && medi->constant=="several statements")
		 //    {
		 //    	medi=medi->children[0];
		 //    }
		 //    r=extract_gamma_coeffs(medi->bracket,1,res,sum);
		 //    if(r!="")
		 //    	res+=r;
			// fout<<sum<<endl;
		 //    fout<<res;
			//fout<<"[#Input current value of pv there.]"<<endl;
    	}
    	else if(constant==":=")
    	{
    		// fout<<"assignment"<<endl;
    		// for(int i=0;i<CFG.size();++i)
    		// {
    		// 	if(label_map[CFG[i].u]==label)
    		// 	{
    		// 		if(label_map[CFG[i].v]!=0)
    		// 			fout<<label_map[CFG[i].v]<<endl;
    		// 		else
    		// 			cout<<"Error: "<<label<<" misses a next label."<<endl;
    		// 	}
    		// }
    		//fout<<extract_coeffs(children[1],1)<<endl;
    		//fout<<label<<endl;
    		extern int assnum;
    		assnum=assnum+1;

    		node *m=children[0];
    		while(m!=NULL and m->type!="pvar")
    			m=m->children[0];
    		node *n=children[1];
    		while(n!=NULL and n->type!="rexpr")
				n=n->children[0];
    		fout<<label<<'@'<<m->constant<<"="<<n->polynom;
    		list<string> svInfo=extract_svariable(n);
    		if(svInfo.empty())
    		{
    			fout<<'@'<<"normal"<<endl;
    		}
    		else
    		{
    			fout<<'@'<<"abnormal";
    			for(list<string>::iterator it=svInfo.begin();it!=svInfo.end();++it)
    			{
    				fout<<'@'<<*it;
    			}
    			//fout<<"\n";
    			for(list<string>::iterator it=svInfo.begin();it!=svInfo.end();++it)
    			{
    				fout<<'@'<<distributions[*it].getType();
    			}
    			//fout<<"\n";
    			for(list<string>::iterator it=svInfo.begin();it!=svInfo.end();++it)
    			{
    				fout<<'@'<<distributions[*it].printInfo();
    			}
    		}
//    		fout<<"[========================================================================]\n";
//    		fout<<"[#Input assignment information(distribution,probability,cost,etc.) there.]\n";
//    		fout<<"[ e.g. abnormal;r;DS;0.25 0.75;1 -1;                                     ]\n";
//    		fout<<"[========================================================================]\n"<<endl;
    		//New add
			//Gamma
			// set<string>::iterator iter=pvars.begin();
		 //    for(;iter!=pvars.end();++iter)
		 //    {
		 //    	gamma_matrix[*iter].clear();
		 //    }
			// int sum=0;
		 //    string res=extract_gamma_coeffs(bracket,1,"",sum);
			// fout<<sum<<endl;
		 //    fout<<res;
    		//fout<<"[#Input current value of pv there.]"<<endl;
    	}
    	else if(constant=="if")
    	{
    		if(children[0]->type=="ndbexpr")
    		{
    			if(children[0]->constant=="prob")
    			{
    				//fout<<"probability"<<endl;
					extern int branchnum;                    
					branchnum=branchnum+1;
		      		string str=children[0]->children[0]->constant;
		      		double i;
		      		i=stod(str);
		      		//fout<<i<<endl;
    				// for(int i=0;i<CFG.size();++i)
		    		// {
		    		// 	if(label_map[CFG[i].u]==label)
		    		// 	{
		    		// 		if(label_map[CFG[i].v]!=0)
		    		// 			fout<<label_map[CFG[i].v]<<" ";
		    		// 	}
		    		// }
		      // 		fout<<endl;
      				//New add
					//Gamma
					// set<string>::iterator iter=pvars.begin();
				 //    for(;iter!=pvars.end();++iter)
				 //    {
				 //    	//gamma_number[*iter]=0;
				 //    	gamma_matrix[*iter].clear();
				 //    }
					// int sum=0;
		   //  		string res=extract_gamma_coeffs(bracket,1,"",sum);
					// fout<<sum<<endl;
				 //    fout<<res;
    				//fout<<"[#Input current value of pv there.]"<<endl;
    			}
    			else if(children[0]->constant=="demon")
    				fout<<"demon"<<endl;
    			else if(children[0]->constant=="angel")
    				fout<<"angel"<<endl;
    			else
    			{
    				fout<<"branching"<<endl;
	    			for(int i=0;i<CFG.size();++i)
			    	{
		    			if(label_map[CFG[i].u]==label)
		    			{
		    				if(label_map[CFG[i].v]!=0)
		    					fout<<label_map[CFG[i].v]<<" ";
		    			}
		    		}
		      		fout<<endl;
		      		//Gamma
					string r="";
					set<string>::iterator iter=pvars.begin();
				    for(;iter!=pvars.end();++iter)
				    {
				    	gamma_matrix[*iter].clear();
				    }
				    string res="";
					int sum = 0;
				    res+=extract_gamma_coeffs(bracket,1,"",sum);
				    r=extract_gamma_coeffs(children[1]->bracket,1,res,sum);
				    if(r!="")
				    	res+=r;
					fout<<sum<<endl;
				    fout<<res;
				    iter=pvars.begin();
				    for(;iter!=pvars.end();++iter)
				    {
				    	gamma_matrix[*iter].clear();
				    }
				    res.clear();
				    sum=0;
			    	res+=extract_gamma_coeffs(bracket,1,"",sum);
			    	r=extract_gamma_coeffs(children[2]->bracket,1,res,sum);
			    	if(r!="")
				    	res+=r;
					fout<<sum<<endl;
				    fout<<res;
    			}
    		}
    		else if(children[0]->type=="star")
    		{
    			fout<<"branching"<<endl;
    			for(int i=0;i<CFG.size();++i)
		    	{
	    			if(label_map[CFG[i].u]==label)
	    			{
	    				if(label_map[CFG[i].v]!=0)
	    					fout<<label_map[CFG[i].v]<<" ";
	    			}
	    		}
	      		fout<<endl;
				//Gamma
				// set<string>::iterator iter=pvars.begin();
			 //    for(;iter!=pvars.end();++iter)
			 //    	gamma_matrix[*iter].clear();
				// int sum=0;
		  //   	string res=extract_gamma_coeffs(bracket,1,"",sum);
				// fout<<sum<<endl;
			 //    fout<<res;
			 //    fout<<sum<<endl;
			 //    fout<<res;
    		}
    	}
    	else if(constant=="tick")
    	{
    		fout<<"tick"<<endl;
    		// for(int i=0;i<CFG.size();++i)
    		// {
    		// 	if(label_map[CFG[i].u]==label)
    		// 	{
    		// 		if(label_map[CFG[i].v]!=0)
    		// 			fout<<label_map[CFG[i].v]<<endl;
    		// 		else
    		// 			cout<<"Error: "<<label<<" misses a next label."<<endl;
    		// 	}
    		// }
   //  		set<string>::iterator iter=pvars.begin();
			// for(;iter!=pvars.end();++iter)
			// {
			// 	gamma_matrix[*iter].clear();
			// }
   //  		//string str=extract_coeffs_tick(children[0],1);
   //  		//fout<<str<<endl;
   //  		node *m=children[0];
   //  		while(m!=NULL and m->type!="rexpr")
			// 	m=m->children[0];
			// fout<<m->polynom<<endl;
   //  		//New add
		 //    //Gamma
			// iter=pvars.begin();
		 //    for(;iter!=pvars.end();++iter)
		 //    {
		 //    	gamma_matrix[*iter].clear();
		 //    }
			// int sum=0;
		 //    string res=extract_gamma_coeffs(bracket,1,"",sum);
			// fout<<sum<<endl;
		 //    fout<<res;
    		//fout<<"[#Input current value of pv there.]"<<endl;
    	}
    	//Added by QXD
    	else if(constant=="skip")
    	{
    		fout<<label<<'@'<<"skip"<<endl;
    		extern int assnum;
    		assnum=assnum+1;
   //  		for(int i=0;i<CFG.size();++i)
   //  		{
   //  			if(label_map[CFG[i].u]==label)
   //  			{
   //  				if(label_map[CFG[i].v]!=0)
   //  					fout<<label_map[CFG[i].v]<<endl;
   //  				else
   //  					cout<<"Error: "<<label<<" misses a next label."<<endl;
   //  			}
   //  		}
   //  		int sum=0;
		 //    string res=extract_gamma_coeffs(bracket,1,"",sum);
			// fout<<sum<<endl;
		 //    fout<<res;
    	}
    }
    // fout<<branchnum<<endl;
    if(type=="expr" or type=="rexpr")
		cout<<"Polynom: |"<<polynom<<"|"<<endl;
	if(type=="polyexpr")
	{
		cout<<"Polynom_set:"<<endl;
		for(set<poly>::iterator iter=polynom_set.begin();iter!=polynom_set.end();++iter)
			cout<<"|"<<*iter<<"|"<<endl;
	}
	if(type=="bexpr")
		set_set_print();
    cout<<"Children: ";
    for(int i=0;i<children.size();++i)
      cout<<children[i]<<"\t";
    if(bracket!=NULL)
      cout<<"Bracket: "<<bracket;
    cout<<endl<<"--------------------"<<endl;
    for(int i=0;i<children.size();++i)
      children[i]->print();
    if(bracket!=NULL)
      bracket->print();
  }

  /*void top_down_pre_eta_calc()
  {
	  if(type!="stmt")
		return;
	  if(constant=="if" and children[0]->constant=="prob")
		{
			int lif=children[1]->label;
			int lelse=children[2]->label;
			string c=children[0]->children[0]->constant;
			stringstream ccnv;
			ccnv<<c;
			double r;
			ccnv>>r;
			//pre_eta[label]=simplify(r*label_polynomial[lif]+(1-r)*label_polynomial[lelse]);
		}
	  else if(constant==":=")
	  {
		  
	  }
	  for(int i=0;i<children.size();++i)
		children[i]->top_down_pre_eta_calc();
  }*/

  //New add
  list<string> extract_svariable(node *expression)
  {
  	list<string> res;
  	if(expression->constant=="single polyexpr" or expression->constant=="single rexpr")
  	{
  		list<string>::iterator it;
  		list<string> c=extract_svariable(expression->children[0]);
		for(it=c.begin();it!=c.end();++it)
			res.push_back(*it);
  	}
  	else if(expression->constant=="+" or expression->constant=="-")
	{
  		list<string>::iterator it;
		list<string> c1=extract_svariable(expression->children[0]);
		for(it=c1.begin();it!=c1.end();++it)
			res.push_back(*it);
		list<string> c2=extract_svariable(expression->children[1]);
		for(it=c2.begin();it!=c2.end();++it)
			res.push_back(*it);
	}
  	else if(expression->constant=="single rvar")
  	{
  		string sv=expression->children[0]->constant;
  		res.push_back(sv);
  	}
  	return res;
  }

  //New add
  string extract_gamma_coeffs(node *expression,int plus_minus,string former_result,int& sum)
  {
  	string res="";
  	string r="";
  	if(expression->constant=="single rexpr")
  	{
		r+=extract_gamma_coeffs(expression->children[0],plus_minus,former_result,sum);
		if(former_result.find(r)==string::npos)
			res+=r;
  	}
  	else if(expression->constant=="single polyexpr")
  	{
		r+=extract_gamma_coeffs(expression->children[0],plus_minus,former_result,sum);
		if(former_result.find(r)==string::npos)
			res+=r;
  	}
  	else if(expression->constant=="and")
  	{
  		set<string>::iterator iter=pvars.begin();
		for(;iter!=pvars.end();++iter)
		{
			gamma_matrix[*iter].clear();
		}
  		string cons=to_string(extract_coeffs(expression->children[0],plus_minus));
  		iter=pvars.begin();
		for(;iter!=pvars.end();++iter)
		{
			if(gamma_matrix[*iter]=="")
			{
				r+="0 ";
			}
			else
			{
				r+=gamma_matrix[*iter]+" ";
			}
		}
		r+=cons;
		if(former_result.find(r)==string::npos)
		{
			res+=(r+"\n");
  			sum++;
		}
		string sub_res=extract_gamma_coeffs(expression->children[1],plus_minus,former_result,sum);
		if(sub_res!="")
  			res+=sub_res;
  	}
  	else
  	{
  		set<string>::iterator iter=pvars.begin();
		for(;iter!=pvars.end();++iter)
		{
			gamma_matrix[*iter].clear();
		}
  		string cons=to_string(extract_coeffs(expression->children[0],plus_minus));
  		iter=pvars.begin();
		for(;iter!=pvars.end();++iter)
		{
			if(gamma_matrix[*iter]=="")
			{
				r+="0 ";
			}
			else
			{
				r+=gamma_matrix[*iter]+" ";
			}
		}
		r+=cons;
		if(former_result.find(r)==string::npos)
		{
			res+=(r+"\n");
  			sum++;
		}
  	}
  	return res;
  }

  double extract_coeffs(node *expression,int plus_minus)
  {
  	double res=0;
  	//cout<<"extract_coeffs: "<<expression->constant<<endl;
  	//cout<<"Type: "<<expression->type<<endl;
  	//cout<<"Result: "<<res<<endl;
  	//fout<<"Extract coefficients"<<endl;
  	if(expression->constant=="single rexpr")
  	{
		res+=extract_coeffs(expression->children[0],plus_minus);
  	}
  	if(expression->constant=="single polyexpr")
  	{
		res+=extract_coeffs(expression->children[0],plus_minus);
  	}
  	if(expression->type=="constant")
	{
		res+=stod(expression->constant)*plus_minus;
		//fout<<stod(child_node->children[0]->constant)*plus_minus<<" ";
	}
	else if(expression->type=="pvar")
	{
		gamma_matrix[expression->constant]=to_string(plus_minus);
	}
	else if(expression->constant=="+")
	{
		res+=extract_coeffs(expression->children[0],plus_minus);
		res+=extract_coeffs(expression->children[1],plus_minus);
	}
	else if(expression->constant=="-")
	{
		res+=extract_coeffs(expression->children[0],plus_minus);
		res+=extract_coeffs(expression->children[1],-plus_minus);
	}
	else if(expression->constant=="single literal" || expression->constant=="()" || expression->constant=="single constant")
	{
		res+=extract_coeffs(expression->children[0],plus_minus);
	}
	else if(expression->constant==">=")
	{
		res+=extract_coeffs(expression->children[0],plus_minus);
	}
	//else if(expression->constant=="<=")
	//{
	//	res+=extract_coeffs(expression->children[1],plus_minus);
	//}
	else if(expression->constant=="*")
	{
		if(expression->children[0]->children[0]->type=="constant" && expression->children[1]->children[0]->type=="pvar")
		{
			gamma_matrix[expression->children[1]->children[0]->constant]=to_string(stod(expression->children[0]->children[0]->constant)*plus_minus);
		}
		else if(expression->children[1]->children[0]->type=="constant" && expression->children[0]->children[0]->type=="pvar")
		{
			gamma_matrix[expression->children[0]->children[0]->constant]=to_string(stod(expression->children[1]->children[0]->constant)*plus_minus);
		}
	}
	else if(expression->constant=="single pvar")
	{
		//cout<<"plus_minus: "<<plus_minus<<endl;
		gamma_matrix[expression->children[0]->constant]=to_string(plus_minus);
	}
	else if(expression->constant=="single rvar")
	{
		//res+=plus_minus;
	}
  	//fout<<"Gamma coefficients:  "<<res<<endl;
  	return res;
  }

};

node *root;

void process_CFG()
{
  cout<<"Process the CFG";
  for(int i=0;i<CFG.size();++i)
    {
	  //cout<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
      edge &e=CFG[i];
      node *u=(node *)e.u,*v=(node *)e.v;
      //cout<<u<<endl;
      //cout<<v<<endl;
      while(u!=NULL and u->constant=="several statements")
	u=u->children[1];
      while(v!=NULL and v->constant=="several statements")
	v=v->children[0];
      e.u = (int *)u;
      e.v = (int *)v;
	  //cout<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
    }
  //process the texts
  for(int i=0;i<CFG.size();++i)
    {
      edge &e=CFG[i];
      if(e.text=="~(demon)")
	e.text="demon";
      if(e.text.substr(0,6)=="~(prob")
	{
		
	   //cout<<'@'<<e.text<<endl;

      int l1;
      int l2;
	  for(int j=0;j<e.text.length();++j)
       {
	    if(isdigit(e.text[j]))
	     {
	    	l1=j;
	        break;
	     }
	   }
	  for(int j=e.text.length()-1;j>-1;--j)
       {
	     if(isdigit(e.text[j]))
	     {
	    	 l2=j;
	    	 break;
	     }
	   }

	      

	  // cout<<'@'<<l1<<endl;
	  // cout<<'@'<<l2<<endl;
	  // cout<<'@'<<e.text<<endl;


	  // stringstream cnv;
	  // cnv<<e.text;
	  // double p;
	  // cnv>>p;

   //    cout<<'3'<<e.text<<endl;
	  // cout<<'4'<<p<<endl;
    
     char ch= '/';

     int loc;
   
     string s1;
     string s2;
     string s3;

     

     //double p;
    
    if(e.text.find(ch)!=string::npos)
    { 
    	//loc1=e.text.find(br1);
    	loc=e.text.find(ch);
    	//loc3=e.text.find(br2);
    	
     // cout<<loc1<<endl;
        //cout<<'$'<<loc<<endl;
      //cout<<loc3<<endl;

         s1=e.text.substr(l1,loc-l1);
         s2=e.text.substr(loc+1,l2-loc);

        // cout<<'@'<<s1<<endl;
        // cout<<'@'<<s2<<endl;

        stringstream ss1;
	    ss1<<s1;
	    double frac1;
	    ss1>>frac1;
	    //cout<<frac1<<endl;

	    stringstream ss2;
	    ss2<<s2;
	    double frac2;
	    ss2>>frac2;
	   // cout<<frac2<<endl;
        double diff;
        diff=frac2-frac1;

        stringstream ss3;
	    ss3<<diff;
	    string numerator;
	    ss3>>numerator;
       // cout<<frac<<endl;
        // double p;
        // p=(1-frac);
        string frac;
        frac=numerator+'/'+s2;

          stringstream cnv2;
	      cnv2<<frac;
	      cnv2>>e.text;
	    // cout<<'#'<<e.text<<'#'<<p<<endl;
	     e.text="prob("+e.text+")";
         }
    else
    {
        //cout<<"NO"<<endl;
         s3=e.text.substr(l1,l2-l1+1);
         //cout<<'@'<<s3<<endl;
         stringstream cnv;
	     cnv<<s3;
	     double p;
	     cnv>>p;

	     p=1-p;
	    // cout<<'p'<<p<<endl;

	     stringstream cnv2;
	  cnv2<<p;
	  cnv2>>e.text;
	 // cout<<'#'<<e.text<<'#'<<p<<endl;
	  e.text="prob("+e.text+")";

    }


   //    cout<<'5'<<p<<endl;

	 

	}
    }
    cout<<"...Done"<<endl;
}

//Processing the information at the head of the program, added by QXD
void process_INFO(int begin,int end)
{
	cout<<"( "<<begin<<" , "<<end<<" )"<<endl;
  	cout<<program<<endl;
  	if(begin<end)
  	{
	  	skip_spaces(begin,end);
	  	int count=0;
	  	int i=begin;
	  	while(i<end)
	  	{
	      if(part(program,i,i+1)=="~")
		  {
		  	int nameEnd=i;
		  	skip_spaces(begin,nameEnd);
    		string variableName=part(program,begin,nameEnd);
    		begin=i+1;
    		int lineEnd;
    		for(int ii=begin;ii<end;ii++)
    			if(program[ii]==';')
    			{
					lineEnd=ii;
					break;
    			}
    		int temp=lineEnd;
	  	    string para="";
	  	    string eval="";
	  	    i=begin;
	  	    if(part(program,i,i+5)=="DUnif")
	    	{
	    		while(begin<lineEnd and program[begin]!='(')
	    			begin++;
	    		while(begin<lineEnd and program[lineEnd]!=')')
	    			lineEnd--;
	  	    	skip_spaces(begin,lineEnd);
	    		for(int ii=begin;ii<lineEnd;ii++)
	    		{
	    			if(program[ii]==',')
	    			{
	    				para=para+part(program,begin+1,ii)+" ";
	    				para=para+part(program,ii+1,lineEnd)+" ";
	    				break;
	    			}
	    		}
	    		distributions[variableName].setType("DU");
	    		distributions[variableName].setPara(para);
	    		distributions[variableName].setEval(eval);
	    	}
	    	else if(part(program,i,i+5)=="CUnif")
	    	{
	    		while(begin<lineEnd and program[begin]!='(')
	    			begin++;
	    		while(begin<lineEnd and program[lineEnd]!=')')
	    			lineEnd--;
	  	    	skip_spaces(begin,lineEnd);
	    		for(int ii=begin;ii<lineEnd;ii++)
	    		{
	    			if(program[ii]==',')
	    			{
	    				para=para+part(program,begin+1,ii)+" ";
	    				para=para+part(program,ii+1,lineEnd)+" ";
	    				break;
	    			}
	    		}
	    		distributions[variableName].setType("CU");
	    		distributions[variableName].setPara(para);
	    		distributions[variableName].setEval(eval);
	    	}
	      	begin=temp+1;
	      	i=begin;
	      	
	  	  }
	  	  else if(part(program,i,i+1)=="P")
	  	  {
    		begin=i+1;
	  	  	string variableName;
	  	  	int lineEnd;
	  	  	for(int ii=begin;ii<end;ii++)
    			if(program[ii]==';')
    			{
					lineEnd=ii;
					break;
    			}
	  	    string para="";
	  	    string eval="";
    		int termEnd;
    		for(int ii=begin;ii<=lineEnd;ii++)
    		{
    			if(program[ii]==';' or program[ii]==',')
    			{
					termEnd=ii;
					int assginEnd=termEnd;
					while(begin<assginEnd and program[begin]!='(')
		    			begin++;
			    	while(begin<assginEnd and program[assginEnd]!=')')
			    		assginEnd--;
			    	skip_spaces(begin,assginEnd);
			  	  	for(int iii=begin;iii<assginEnd;iii++)
		    		{
		    			if(program[iii]=='=')
		    			{
		    				variableName=part(program,begin+1,iii);
		    				if(distributions.find(variableName)!=distributions.end())
		    				{
		    					distributions.find(variableName)->second.addPara(part(program,iii+1,assginEnd));
		    				}
		    				else
		    				{
		    					para=para+part(program,iii+1,assginEnd)+" ";
	    						distributions[variableName].setType("DS");
					    		distributions[variableName].setPara(para);
					    		distributions[variableName].setEval(eval);
		    				}
		    				break;
		    			}
		    		}
		    		begin=assginEnd+1;
		    		int probEnd=termEnd;
		    		skip_spaces(begin,probEnd);
		    		for(int iii=begin;iii<probEnd;iii++)
		    		{
		    			if(program[iii]=='=')
		    			{
		    				begin=iii+1;
		    				if(distributions.find(variableName)!=distributions.end())
		    					distributions.find(variableName)->second.addEval(part(program,begin,probEnd));
		    				else
		    					cout<<"Error: "<<variableName<<" should already be stored."<<endl;
		    				break;
		    			}
		    		}
	      			begin=termEnd+1;
	      			ii=begin;
    			}
    		}
    		begin=lineEnd;
    		i=begin;
	  	  }
	  	  i++;
	  	}
	  	fout<<"linear"<<endl;
  	}
  	else
  	{   
  		string stype;
  		stype="affine";
        fout<<stype<<endl;
  		cout<<"No sampling variables."<<endl;
  		//cout<<"Error: Infomation can not be parsed. Begin is "<<begin<<". End is "<<end<<"."<<endl;
  	}
}

int main(int argc, char* argv[])
{
 extern int branchnum;
 extern int assnum;
  if(argc<3)
  {
  	cout<<"Not enough parameters."<<endl;
  }
  else if(argc>3)
  {
  	cout<<"Too much parameters."<<endl;
  }
  else
  {
  	fout.open(argv[1]);
  	cout<<"Output file 1: "<<argv[1]<<endl;
  	cout<<"Output file 2: "<<argv[2]<<endl;
  }
  int r,i;

  //New Add
  //Parse the distribution of sampling variables by QXD
  for(i=0;(r=getchar())!='#';i++)
  	input[i]=r;
  program=input;
  memset(input,0,sizeof(input));
  cout<<"Parse Infomation:"<<endl;
  process_INFO(0,program.length());
  unordered_map<string,info>::iterator it=distributions.begin();
  for(;it!=distributions.end();++it)
  {
  	cout<<"Sampling variable: "<<it->first<<endl;
  	cout<<it->second.printInfo()<<endl;
  }
  //New Add
  //Clear the "#" clause at the beginning of each inputs
   //cout<<"Sensitivity Type : "<<endl;
  // for(i=0;(r=getchar())!='#';i++)
  // {cout<<char(r);fout<<char(r);}
  // cout<<"\n";

  cout<<"Loop guards: "<<endl;
  for(i=0;(r=getchar())!='#';i++)
  {
  	cout<<char(r);
  	fout<<char(r);

 //  	set<string>::iterator iter=pvars.begin();
 //    for(;iter!=pvars.end();++iter)
 //    {
 //    	//gamma_number[*iter]=0;
 //    	gamma_matrix[*iter].clear();
 //    }
	// int sum=0;
	// string res=extract_gamma_coeffs(bracket,1,"",sum);
	// fout<<sum<<endl;
 //    fout<<res;
}
  cout<<"\n";

   cout<<"Violated guards: "<<endl;
  for(i=0;(r=getchar())!='#';i++)
  {cout<<char(r);fout<<char(r);}
  cout<<"\n";
  //========
  for(i=0;(r=getchar())!=EOF;i++)
    input[i]=r;
  input[i]=0;
  program=input;
  cout<<"Input Code:"<<endl;
  cout<<program<<endl;
  cout<<"Parse Tree:"<<endl;
  root=new node("stmt",0,program.length());
  fout<<pvars.size()<<endl;
  set<string>::iterator iter=pvars.begin();
  for(;iter!=pvars.end();++iter)
  {
  	fout<<*iter<<" ";
  }
  fout<<endl;
  //fout<<"[#Input degree d there.]"<<endl;
  fout<<1<<endl;
  //fout<<"[#Input k there.]"<<endl;
  fout<<1<<endl;
 
  for(int i=0;i<CFG.size();++i)
  {
  	if(label_map[CFG[i].u]>max_label)
  		max_label=label_map[CFG[i].u];
  	if(label_map[CFG[i].v]>max_label)
  		max_label=label_map[CFG[i].v];
  }
  //fout<<max_label<<endl;
  //fout<<"[#Input matrix of program variables value there.]"<<endl;
  //fout<<"[#Input ininitial value of program variables value there.]"<<endl;
  process_CFG();
  root->print();
  cout<<endl<<"CFG:"<<endl;
  for(int i=0;i<CFG.size();++i)
    cout<<label_map[CFG[i].u]<<" -> "<<label_map[CFG[i].v]<<" :: "<<CFG[i].text<<endl;
  fout<<"out"<<endl;
  fout.close();

fout.open(argv[2]);
  //fout<<branchnum<<endl;
fout<<max_label<<endl;
fout<<CFG.size()<<endl;
 for(int i=0;i<CFG.size();++i)
    fout<<label_map[CFG[i].u]<<"*"<<label_map[CFG[i].v]<<"*"<<CFG[i].text<<endl;

// for(int i=0;i<CFG.size();++i)
//     fout<<CFG[i].u<<" -> "<<CFG[i].v<<endl;
  //create the monomials
  //root->recursively_create_monomials();
  
  //calculate pre_etas
  //root->recursively_calculate_pre_etas();  
  
  //write the final output
  //root->recursively_write_output();
  return 0;
}
