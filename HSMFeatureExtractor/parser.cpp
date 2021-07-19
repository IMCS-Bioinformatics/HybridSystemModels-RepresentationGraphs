#include <bits/stdc++.h>

using namespace std;

// for each gene its function in string form
map<string, string> functions;
map<string, vector<string>> f_arguments;

// list of genes
vector<string> genes;

vector<string> bs_succession;
// for each bindingsite what could be attached 
map<string,vector<string>> bindingsiteregulators;

// this holds for each binding site what is attached (empty string if nothing)
map<string, string> current_state;

// constraints
vector<string> constraints;

bool variable_value(string var) {
	// gets variable, looks at which state it is at the moment and
	// returns 1 if it is present at the current state and 0 otherwise
	int pos1 = var.find("[");
	int pos2 = var.find("]");
	string checkpos =  var.substr(pos1+1, pos2 - pos1-1);
	string variable = var.substr(0, pos1);
	bool ans = (current_state[variable] == checkpos);

	return ans;
}


// given a function calculates its value in given current state
// assumes basic boolean algebra, + is disjunction, * is conjunction and ~ is negation
bool calculate(string clause) {

	int case_difference = 0;
	int position = -1;
	// does it containt disjunction outside of parentheses?
	for (int i = 0; i < clause.size(); i++) {
		if (clause[i] == '(') {
			case_difference++;
		}
		else if (clause[i] == ')') {
			case_difference--;
		}
		else if (clause[i] == '+' && case_difference == 0) {
			position = i;
			break;
		}
	}

	if (position > 0) {
		return calculate(clause.substr(0,position)) | calculate(clause.substr(position+1));
	}

	case_difference = 0;
	position = -1;
	// does it containt conjunction outside of parentheses?
	for (int i = 0; i < clause.size(); i++) {
		if (clause[i] == '(') {
			case_difference++;
		}
		else if (clause[i] == ')') {
			case_difference--;
		}
		else if (clause[i] == '*' && case_difference == 0) {
			position = i;
			break;
		}
	}

	if (position > 0) {
		return calculate(clause.substr(0,position)) & calculate(clause.substr(position+1));
	}


	// if gets this far - either is ~(), () or a variable
	// 
	if (clause[0] == '~') {
		return !calculate(clause.substr(1));
	}

	if (clause[0] == '(' && clause[clause.size()-1] == ')') {
		return calculate(clause.substr(1,clause.size()-2));
	}

	
	return variable_value(clause);
}




int main(int argc, char **argv) {
	// reading in
	if (argc < 3) {
		printf("Need filenames\n");
		return 0;
	}



	//reading in the model

	ifstream modelInfo;
	string filename = (string) argv[1];
	modelInfo.open(filename);
	string line;

	ofstream output;
	string outfilename = (string) argv[2];



	int current_case = 0;
	// 1 - genes
	// 2 - binding sites
	// 3 - regulations and functions
	// 4 - constraints_affinities


	while (getline(modelInfo,line)){		
		// ignore comments
		if (line[0] == '%') {
			continue;
		}

		// ignore whitespace lines
		size_t pos = line.find_first_not_of(" 	");
		if (pos >= line.size()) {
			continue;
		}
		
		line = line.substr(pos);
		// which part of the model are we looking at?
		if (line[0] == '#') {
			if (line.find("#GENES") != string::npos) {
				current_case = 1;
			}
			else if (line.find("#REGULATIONS") != string::npos) {
				current_case = 3;
			}
			else if (line.find( "#BINDINGSITES") != string::npos) {
				current_case = 2;
			}
			else if (line.find("#CONSTRAINTS_AFFINITIES") != string::npos){
				current_case = 4;
			}
			else if (line[0] == '#'){
				current_case = -1;
				cout << line << endl;
			}
			continue;
		}
	

		string gene;
		string bindingsite;
		while(line[line.size()-1] != ';' && line.size() > 0) {
			line.pop_back();
		}

		if (line.size() == 0) {
			continue;
		}

		switch(current_case) {
			case 1:
				//	split line on spaces
				//	add genes to the list 
				if (line[line.size()-1] != ';') {
					cout << "Lines have to end in ;\n";
					cout << line;
					return 0;
				}

				line.pop_back();
				gene = "";

				for (auto i: line) {
					if (i != ' ' & i != '	') {
						gene += i;
					}
					else if (gene != "") {
						genes.push_back(gene);
						gene = "";
					}
				}
				genes.push_back(gene);
				break;

			case 2:
				if (line[line.size()-1] != ';') {
					cout << "Lines have to end in ;\n";
					return 0;
				}

				line.pop_back();
				bindingsite = "";

				for (auto i: line) {
					if (i != ' ' & i != '	') {
						bindingsite += i;
					}
					else if (bindingsite != "") {

						int pos1 = bindingsite.find('{');
						int pos2 = bindingsite.find('}');
						vector<string> sites;
						sites.push_back("");
						string tmp = "";
						for (auto b: bindingsite.substr(pos1+1, pos2 - pos1 - 1)) {
							if (b != ',') {
								tmp += b;
							}
							else {
								sites.push_back(tmp);
								tmp = "";
							}
						}
						sites.push_back(tmp);
						string bs = bindingsite.substr(0, pos1); 
						bindingsiteregulators[bs] = sites;
						bs_succession.push_back(bs);
						bindingsite = "";
					}
				}
 					if (bindingsite != "") {

						int pos1 = bindingsite.find('{');
						int pos2 = bindingsite.find('}');
						vector<string> sites;
						sites.push_back("");
						string tmp = "";
						for (auto b: bindingsite.substr(pos1+1, pos2 - pos1 - 1)) {
							if (b != ',') {
								tmp += b;
							}
							else {
								sites.push_back(tmp);
								tmp = "";
							}
						}
						sites.push_back(tmp);
						string bs = bindingsite.substr(0, pos1); 
						bindingsiteregulators[bs] = sites;
						bs_succession.push_back(bs);
						bindingsite = "";
					}
					
					
				break;

			case 3:
				// read and add in functions
				

				string line2 = "";
				// remove whitespaces

				for (auto c: line) {
					if (c != ' ' && c != '	') {
						line2 += c;
					}
				}

				
				//split on '='
				string first = line2.substr(0, line2.find("="));
				string second = line2.substr(line2.find("=")+1);
				second = second.substr(0, second.size()-1);

				// everything in cases '(',')' before '='  seperated by commas are variables 
				vector<string> variables;                                                                                    
				string name = first.substr(0, first.find("("));
				string varlist = first.substr(first.find("(")+1, first.size() - first.find("(") - 2);

				
				int cpos = 0;
				int ccount = 0;
				for (int i = 0; i < varlist.size(); i++) {
					if (varlist[i] == ',') {
						variables.push_back(varlist.substr(cpos, ccount));
						cpos = i + 1;
						ccount = 0;
					}
					else
						ccount++;
				}
				variables.push_back(varlist.substr(cpos));

				//everything before cases are gene names, seperated by commas
				// if no commas present, only single gene

				if (name.find(",") != -1) {
					while (name.find(",") != string::npos) {
						int pos = name.find(",");
						functions[name.substr(0, pos)] = second;
						f_arguments[name.substr(0, pos)] = variables;
						name = name.substr(pos+1);
					}
				}
				functions[name] = second;
				f_arguments[name] = variables;
				break;

			case 4:
				// read and add in constraints affinities

				constraints.push_back(line);
				break;

			default: 
				cout << "Line '"<<line <<  "' does not fit the supported format." << endl;
		}
	
	}


	cout << "Done parsing!" << endl;

	output.open(outfilename);

	/**/
	//Modelname
	output << filename.substr(0, filename.find_last_of("."))  << endl;

	// output number of genes
	output << genes.size() << endl;

	// calculate number of transcript factors - genes present as arguments and
	// count binary and ternary BS


	vector<string> genes2;

	int countbin = 0, counttern = 0;
	set<string> transcript_factors;
	
	for ( auto const& i : bindingsiteregulators) {
		if ( i.second.size() == 2)
			countbin++;
		else counttern++;
		for (auto tf : i.second) {
			if(tf != "")
				transcript_factors.insert(tf);
		}
	}

	for (auto tf: transcript_factors) {
		genes2.push_back(tf);
	}

	// output genes
	for (auto g: genes) {
		if (transcript_factors.find(g) == transcript_factors.end()) {
			genes2.push_back(g);
		}
	}


	map<string, int> genenr;
	map<int, string> nrgene;

	int nr = 0;
	for (auto g: genes) {
		genenr[g] = nr;
		nrgene[nr] = g;
		nr++;
		output << g << " ";
	}
	output << endl;


	// output number of transcipt factors
	output << transcript_factors.size() << endl;
	
	// output count of binary BS
	output << countbin << endl;;

	// output binary BS
	map<pair<string,string>, int> bindingsitenr;
	map<int, string> nrbindingsite;

	nr = 0;
	for (int j = 0; j < bs_succession.size(); j++) {
		auto i = bindingsiteregulators[bs_succession[j]];
		if (i.size() == 2) {
			output << bs_succession[j] << "{" << i[1] << "} ";
			pair<string,string> pp = make_pair(bs_succession[j], i[1]); 
			bindingsitenr[pp] = nr;
			nrbindingsite[nr] = bs_succession[j];
			nr++;
		}
	}

	output << endl;

	//output count of ternary bs
	output << counttern << endl;
	//output ternary bs
	for (int j = 0; j < bs_succession.size(); j++) {
		auto i = bindingsiteregulators[bs_succession[j]];
		if (i.size() != 2){
			output << bs_succession[j] << "{" << i[1] <<  ","<< i[2]<< "} ";
			bindingsitenr[make_pair(bs_succession[j],i[1])] = nr;
			nrbindingsite[nr] = bs_succession[j];
			nr++;
			bindingsitenr[make_pair(bs_succession[j],i[2])] = nr;
			nrbindingsite[nr] = bs_succession[j];
			nr++;

		}
	}

	output << endl;

	// mapping BS -> genes: 1 Bs for each binary site (first) (cII), 2 for each ternary site (3x (cI & cro)

	for (int i = 0; i < nrbindingsite.size(); ) {

		string bs = nrbindingsite[i];
		vector<string> regulators = bindingsiteregulators[bs];
		for (auto r : regulators) {
			if (r != "")
			output << genenr[r] << " ", i++;
		}


	}
	output << endl;
	

	// # of args for gene reg functions (order of genes: cI,CII,...), list length = # of gene

	for (string g: genes) {
		int tmp = 0;
		for (auto i : f_arguments[g]) {
			tmp += bindingsiteregulators[i].size() - 1;
			for (auto j : bindingsiteregulators[i]) {
			}
		}
		output <<  tmp << " ";

	}
	output << endl << endl;


	// for each gene its binding sites
	for (string g: genes) {


		for (int i = 0; i <  f_arguments[g].size(); i++) {
			
			for (int p =   1; p < bindingsiteregulators[f_arguments[g][i]].size(); p++)
			
			output << 
				bindingsitenr[
					make_pair(
						f_arguments[g][i],
						bindingsiteregulators[f_arguments[g][i]][p]
					)]
				<< " ";
		}
		output << endl;
	}


    // for each gene calculate its function's values               
	for (string g: genes) {
		string result = "";
		
		//calculate length of result
		
		int dimensions = 0;
		
		for (int i = 0; i < f_arguments[g].size(); i++) {
			dimensions += bindingsiteregulators[f_arguments[g][i]].size() -1;
		}
		


		for (long long i = 0; i < 1 << dimensions; i++){
			
			//prepare state depending on bitmask

			long long iterator = i;
			int position = 0;
			int position_vec = 0;
			for (int j = f_arguments[g].size()-1; j >= 0; j--) {
				current_state[f_arguments[g][j]] = "";
				
				for (int p = bindingsiteregulators[f_arguments[g][j]].size()-1; p >0; p--) {

				
					if ( iterator % 2 == 1) {
						current_state[f_arguments[g][j]] = bindingsiteregulators[f_arguments[g][j]][p];
					}
					iterator /= 2;
				}
			}

			//calculate the value of function in state
			//and append this to string of function values
			bool tmp = calculate(functions[g]);

			if (tmp == false) {
				result += "0 ";
			}
			else 
				result += "1 ";
		}

		output << result << endl;
	}

	/*
	// use this to output the read in constraints
	// output constraints
	output << constraints.size() << endl;
	for (auto c: constraints) {
		output << c << endl;
	}
	*/
	/**/
	cout << "Done outputting" << endl;

	return 0;
	
}