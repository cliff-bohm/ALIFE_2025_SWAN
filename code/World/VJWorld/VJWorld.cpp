//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "VJWorld.h"

// this is how you setup a parameter in MABE, the function Parameters::register_parameter()takes the
// name of the parameter (catagory-name), default value (which must conform with the type), a the useage message
shared_ptr<ParameterLink<int>> VJWorld::evaluationsPerGenerationPL =
    Parameters::register_parameter("WORLD_VJ-evaluationsPerGeneration", 30,
    "how many times should each organism be tested in each generation?");

// the constructor gets called once when MABE starts up. use this to set things up
VJWorld::VJWorld(shared_ptr<ParametersTable> PT) : AbstractWorld(PT) {
    
    //localize a parameter value for faster access
    evaluationsPerGeneration = evaluationsPerGenerationPL->get(PT);
    
    std::fill(more_ones.begin(), more_ones.begin() + 55, 1);
    std::fill(more_ones.begin() + 55, more_ones.end(), 0);
    std::fill(less_ones.begin(), less_ones.begin() + 45, 1);
    std::fill(less_ones.begin() + 45, less_ones.end(), 0);

    // popFileColumns tell MABE what data should be saved to pop.csv files
	popFileColumns.clear();
    popFileColumns.push_back("score");
    popFileColumns.push_back("less_score");
    popFileColumns.push_back("more_score");
}

// the evaluate function gets called every generation. evaluate should set values on organisms datamaps
// that will be used by other parts of MABE for things like reproduction and archiving
auto VJWorld::evaluate(map<string, shared_ptr<Group>>& groups, int analyze, int visualize, int debug) -> void {

    int o0 = 0;
    int o1 = 0;
    bool done = false;
    int popSize = groups[groupName]->population.size(); 
    
    // evaluate this organism some number of times based on evaluationsPerGeneration
    for (int eval = 0; eval < evaluationsPerGeneration; eval++) {
        //std::cout << eval << "   eval number" << std::endl;
        std::shuffle(more_ones.begin(), more_ones.end(), Random::getCommonGenerator());
        std::shuffle(less_ones.begin(), less_ones.end(), Random::getCommonGenerator());

        // in this world, organisms donot interact, so we can just iterate over the population
        for (int i = 0; i < popSize; i++) {
            // create a shortcut to access the organism and organisms brain
            auto org = groups[groupName]->population[i];
            auto brain = org->brains[brainName];
            brain->resetBrain();
            done = false;
            for (int j = 0; j < 100 && !done; j++) {
                brain->setInput(0, more_ones[j]);
                brain->update();
                //std::cout << " " << j;
                if (j >= 80) {
                    o0 = Bit(brain->readOutput(0));
                    o1 = Bit(brain->readOutput(1));
                    //std::cout << "* " << o0 << o1;
                    if (o0 > o1) {
                        //std::cout << "M";
                        org->dataMap.append("score", 1.0);
                        org->dataMap.append("more_score", 1.0);
                        done = true;
                    }
                    if (o1 > o0) {
                        //std::cout << "l";
                        org->dataMap.append("score", 0.0);
                        org->dataMap.append("more_score", 0.0);
                        done = true;
                    }
                }
            }
            if (!done) {
                org->dataMap.append("score", 0.0);
                org->dataMap.append("more_score", 0.0);
            }
            brain->resetBrain();
            done = false;
            for (int j = 0; j < 100 && !done; j++) {
                brain->setInput(0, less_ones[j]);
                brain->update();
                //std::cout << " " << j;
                if (j >= 80) {
                    o0 = Bit(brain->readOutput(0));
                    o1 = Bit(brain->readOutput(1));
                    //std::cout << "* " << o0 << o1;
                    if (o0 > o1) {
                        //std::cout << "m";
                        org->dataMap.append("score", 0.0);
                        org->dataMap.append("less_score", 0.0);
                        done = true;
                    }
                    if (o1 > o0) {
                        //std::cout << "l";
                        org->dataMap.append("score", 1.0);
                        org->dataMap.append("less_score", 1.0);
                        done = true;
                    }
                }
            }
            if (!done) {
                //std::cout << "NO";
                org->dataMap.append("score", 0.0);
                org->dataMap.append("less_score", 0.0);
            }
            //std::cout << std::endl;
        }
    } // end of population loop
}

// the requiredGroups function lets MABE know how to set up populations of organisms that this world needs
auto VJWorld::requiredGroups() -> unordered_map<string,unordered_set<string>> {
	return { { groupName, { "B:"+brainName+",1,2" } } };
        
        // this tells MABE to make a group called "root::" with a brain called "root::" that takes 2 inputs and has 1 output
        // "root::" here also indicates the namespace for the parameters used to define these elements.
        // "root::" is the default namespace, so parameters defined without a namespace are "root::" parameters
}
