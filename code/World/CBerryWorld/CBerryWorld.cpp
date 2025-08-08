//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "CBerryWorld.h"

long long CBerryWorld::Plant::nextID = 0;
long long CBerryWorld::Critter::nextID = 0;

// this is how you setup a parameter in MABE, the function Parameters::register_parameter()takes the
// name of the parameter (catagory-name), default value (which must conform with the type), a the useage message
std::shared_ptr<ParameterLink<int>> CBerryWorld::evaluationsPerGenerationPL =
    Parameters::register_parameter("WORLD_CBerry-evaluationsPerGeneration", 30,
    "how many times should each organism be tested in each generation?");




// the constructor gets called once when MABE starts up. use this to set things up
CBerryWorld::CBerryWorld(std::shared_ptr<ParametersTable> PT) : AbstractWorld(PT) {

    //localize a parameter value for faster access
    evaluationsPerGeneration = evaluationsPerGenerationPL->get(PT);
    
    // popFileColumns tell MABE what data should be saved to pop.csv files
	popFileColumns.clear();
    popFileColumns.push_back("score");
    popFileColumns.push_back("moves");
    popFileColumns.push_back("turns");
    popFileColumns.push_back("eats");
}

// the evaluate function gets called every generation. evaluate should set values on organisms datamaps
// that will be used by other parts of MABE for things like reproduction and archiving
auto CBerryWorld::evaluate(std::map<std::string, std::shared_ptr<Group>>& groups, int analyze, int visualize, int debug) -> void {

    //
    // set up world
    // 
    CGrid<Sector> world;
    std::vector<std::shared_ptr<Plant>> plants;
    std::vector< std::shared_ptr<Critter>> critters;

    std::vector<std::array<CPoint, 4>> viewWindows;
    viewWindows.push_back(makeViewWindow(.2, .4, 10, 0));
    viewWindows.push_back(makeViewWindow(.4, .7, 10, 0));
    viewWindows.push_back(makeViewWindow(.7, 1.0, 30, 0));
    viewWindows.push_back(makeViewWindow(.2, .7, 20, 15));
    viewWindows.push_back(makeViewWindow(.2, .7, 20, -15));
    //viewWindows.push_back(makeViewWindow(.2, .7, 40, 35));
    //viewWindows.push_back(makeViewWindow(.2, .7, 40, -35));
    //viewWindows.push_back(makeViewWindow(.7, 1, 40, 35));
    //viewWindows.push_back(makeViewWindow(.7, 1, 40, -35));

    if (Global::update == 0) {
        int c = 0;
        for (auto& vw : viewWindows) {
            std::cout << "x_" << std::to_string(c) << " = [";
            for (auto& p : vw) {
                std::cout << p.x << ",";
            }
            std::cout << "]" << std::endl << "y_" << std::to_string(c) << " = [";
            for (auto& p : vw) {
                std::cout << p.y << ",";
            }
            std::cout << "]" << std::endl;
            c += 1;
        }
    }

    world.reset(worldX, worldY);
    for (int x = 0; x < worldX; x++) {
        for (int y = 0; y < worldY; y++) {
            world(x, y).loc = CPoint(x, y);
        }
    }

    // place the plants
    int x = 0;
    int y = 0;

    for (auto& plant : groups["plant::"]->population) {
        plants.push_back(std::make_shared<Plant>());
        plants.back()->org = plant;
        plants.back()->loc = CPoint(x + Random::getDouble(1.0),y + Random::getDouble(1.0));
        world(int(plants.back()->loc.x), int(plants.back()->loc.y)).plants.push_back(plants.back());
        plants.back()->localIndex = world(plants.back()->loc).plants.size() - 1;
        x += 3;
        if (x >= worldX) {
            y += 1;
            x = loopMod(x,worldX);
            if (y >= worldY) {
                y = 0;
            }
        }
    }

    // place the critters
    for (auto& critter : groups["critter::"]->population) {
        critters.push_back(std::make_shared<Critter>());
        critters.back()->org = critter;
        critters.back()->loc = CPoint(Random::getDouble(worldX), Random::getDouble(worldY));
        //critters.back()->facing = Random::getDouble(360);
        world(critters.back()->loc).critters.push_back(critters.back());
        critters.back()->localIndex = world(critters.back()->loc).critters.size() - 1;
    }

    if (0) {// report
        int plantCount = 0;
        int critterCount = 0;

        for (int x = 0; x < worldX; x++) {
            for (int y = 0; y < worldY; y++) {
                //std::cout << "in sector: " << world(x, y).loc.as_string() << std::endl;
                //std::cout << "  plants count: " << std::to_string(world(x, y).plants.size()) << std::endl;
                for (auto& plant : world(x, y).plants) {
                    //std::cout << "    plant: " << plant->ID << "  " << plant->loc.as_string() << std::endl;
                    plantCount++;
                }
                //std::cout << "  critters count: " << std::to_string(world(x, y).critters.size()) << std::endl;
                for (auto& critter : world(x, y).critters) {
                    //std::cout << "    critter: " << critter->ID << "  " << critter->loc.as_string() << std::endl;
                    critterCount++;
                }
            }
        }
        std::cout << "plantCount: " << plantCount << "   plants.size(): " << plants.size() << std::endl;
        std::cout << "critterCount: " << critterCount << "   critters.size(): " << critters.size() << std::endl;
    }

    for (auto& critter : critters) {


        for (auto& plant : plants) {
            if (!plant->alive) {
                plant->alive = true;
                plant->leafCoverage = 1;
            }
        }


        for (int u = 0; u < 200; u++) { // world updates
            //std::shuffle(std::begin(critters), std::end(critters), Random::Generator());


            std::vector<std::shared_ptr<Critter>> killList; // orgs that die this update
            //for (auto& critter : critters) {

                // rotate each arc to critters facing
            std::vector<std::array<CPoint, 4>> actual_windows(viewWindows.size());
            for (int i = 0; i < viewWindows.size(); i++) {
                actual_windows[i] = rotateViewWindow(viewWindows[i], critter->facing, critter->loc);
            }

            if (0) {
                if (Global::update == 0) {
                    int c = 0;
                    for (auto& vw : actual_windows) {
                        std::cout << "x_" << std::to_string(c) << " = [";
                        for (auto& p : vw) {
                            std::cout << p.x << ",";
                        }
                        std::cout << "]" << std::endl << "y_" << std::to_string(c) << " = [";
                        for (auto& p : vw) {
                            std::cout << p.y << ",";
                        }
                        std::cout << "]" << std::endl;
                        c += 1;
                    }
                }
            }

            double leavesHere = 0;
            double crittersHere = 0;
            std::vector<std::shared_ptr<Plant>> plantsHere;
            std::vector<double> leafTotals(actual_windows.size(), 0);
            std::vector<double> crittersTotals(actual_windows.size(), 0);

            int rx, ry;
            for (int x = -1; x < 2; x++) {
                for (int y = -1; y < 2; y++) {
                    rx = loopMod(critter->loc.x + x, worldX);
                    ry = loopMod(critter->loc.y + y, worldY);
                    for (auto& plant : world(rx, ry).plants) {
                        if (plant->loc.InCircle(critter->loc, .2)) {
                            leavesHere += plant->leafCoverage;
                            plantsHere.push_back(plant);
                            //std::cout << "*";
                        }
                        for (int i = 0; i < actual_windows.size(); i++) {
                            //std::cout << ":";
                            if (plant->loc.InConvex4Poly(actual_windows[i])) {
                                leafTotals[i] += plant->leafCoverage;
                                //std::cout << "+";
                            }
                        }
                        //std::cout << std::endl;
                    }
                    //for (auto& critter : world(rx, ry).critters) {
                    //    if (critter->loc.InCircle(critter->loc, .2)) {
                    //        crittersHere += 1;
                    //    }
                    //    for (int i = 0; i < actual_windows.size(); i++) {
                    //        if (critter->loc.InConvex4Poly(actual_windows[i])) {
                    //            crittersTotals[i] += 1;
                    //        }
                    //    }
                    //}
                }
            }


            int inputCounter = 0;
            critter->org->brains["critter::"]->setInput(inputCounter++, leavesHere);
            for (int i = 0; i < viewWindows.size(); i++) {
                critter->org->brains["critter::"]->setInput(inputCounter++, leafTotals[i]);
            }


            critter->org->brains["critter::"]->update();


            if (critter->org->brains["critter::"]->readOutput(0) > 0) { // eat
                critter->energy += leavesHere;
                //std::cout << "eat: " << critter->energy << "  " << leavesHere << std::endl;
                for (auto& plant : plantsHere) {
                    critter->eatCnt++;
                    plant->alive = false; // he killed it! will clean up later
                    plant->leafCoverage = 0;
                }
            }
            else { // move

                //std::cout << "  " << critter->facing << "  " << critter->loc.as_string() << "  ";
                double left = critter->org->brains["critter::"]->readOutput(1);
                double right = critter->org->brains["critter::"]->readOutput(2);
                double turn = (std::max(0.0, std::min(1.0, right)) - std::max(0.0, std::min(1.0, left))) * 45.0;
                double speed = (std::max(0.0, std::min(1.0, critter->org->brains["critter::"]->readOutput(3))));
                if (speed > 0) {
                    critter->moveCnt++;
                }
                if (turn < -10 || turn > 10) {
                    critter->turnCnt++;
                }
                critter->facing = loopModDouble(critter->facing + turn, 360.0);
                if (world(critter->loc).critters.size() > 1) {
                    world(critter->loc).critters.back()->localIndex = critter->localIndex;
                    std::swap(world(critter->loc).critters[critter->localIndex], world(critter->loc).critters.back());
                }
                world(critter->loc).critters.pop_back();
                auto movement = CPoint(0, speed).rotate(critter->facing);
                critter->loc.x = loopModDouble(critter->loc.x + movement.x, (double)worldX);
                critter->loc.y = loopModDouble(critter->loc.y + movement.y, (double)worldY);
                world(critter->loc).critters.push_back(critter);
                critter->localIndex = world(critter->loc).critters.size() - 1;
                //std::cout << " -> " << critter->facing << "  " << critter->loc.as_string() << "  ";

                //std::cout << "r:" << right << "  l:" << left << "  s:" << speed << std::endl;
            }
        }

        // remove dead critters

        for (auto& critter : critters) { // make babies
            // if this critter is reproducting and has enough energy, make an offspring
        }

        if (0) {
            for (auto& plant : plants) {
                if (!plant->alive) {
                    if (world(plant->loc).plants.size() > 1) {
                        world(plant->loc).plants.back()->localIndex = plant->localIndex;
                        std::swap(world(plant->loc).plants[plant->localIndex], world(plant->loc).plants.back());
                    }
                    world(plant->loc).plants.pop_back();
                    plant->loc = CPoint(Random::getDouble(worldX), Random::getDouble(worldY));
                    world(plant->loc).plants.push_back(plant);
                    plant->localIndex = world(plant->loc).plants.size() - 1;
                    plant->leafCoverage = 1;
                }
                // (later) collect resource and water
                // (later) set plant inputs
                // (later) brain->update
                // (later) update state / grow
            }
        }
        // for now, if plant has no leaves, kill it

        // for now, just make some new plants!
    }// end of population loop

    if (visualize) {
        std::cout << "  score: " << critters[0]->energy  << " + " << .1 * std::min(critters[0]->moveCnt, std::min(critters[0]->turnCnt, critters[0]->eatCnt)) << std::endl;
    }

    for (int i = 0; i < critters.size(); i++) {
        groups["critter::"]->population[i]->dataMap.append("score", critters[i]->energy + .1 * std::min(critters[i]->moveCnt, std::min(critters[i]->turnCnt, critters[i]->eatCnt)));
        groups["critter::"]->population[i]->dataMap.append("moves", critters[i]->moveCnt + .1 * std::min(critters[i]->moveCnt, std::min(critters[i]->turnCnt, critters[i]->eatCnt)));
        groups["critter::"]->population[i]->dataMap.append("turns", critters[i]->turnCnt + .1 * std::min(critters[i]->moveCnt, std::min(critters[i]->turnCnt, critters[i]->eatCnt)));
        groups["critter::"]->population[i]->dataMap.append("eats", critters[i]->eatCnt + .1 * std::min(critters[i]->moveCnt, std::min(critters[i]->turnCnt, critters[i]->eatCnt)));
    }
    for (int i = 0; i < plants.size(); i++) {
        groups["plant::"]->population[i]->dataMap.append("score", plants[i]->energy);
        groups["plant::"]->population[i]->dataMap.append("moves", 0);
        groups["plant::"]->population[i]->dataMap.append("turns", 0);
        groups["plant::"]->population[i]->dataMap.append("eats", 0);
    }

    critters.clear();
    plants.clear();

    for (int x = 0; x < worldX; x++) {
        for (int y = 0; y < worldY; y++) {
            world(x, y).plants.clear();
            world(x, y).critters.clear();
        }
    }
}

// the requiredGroups function lets MABE know how to set up populations of organisms that this world needs
auto CBerryWorld::requiredGroups() -> std::unordered_map<std::string, std::unordered_set<std::string>> {
    return { { (std::string)"plant::", {"B:" + (std::string)"plant::" + ",2,1"}},
        {(std::string)"critter::", {"B:" + (std::string)"critter::" + ",6,4"}} };
        
        // this tells MABE to make a group called "root::" with a brain called "root::" that takes 2 inputs and has 1 output
        // "root::" here also indicates the namespace for the parameters used to define these elements.
        // "root::" is the default namespace, so parameters defined without a namespace are "root::" parameters
}
