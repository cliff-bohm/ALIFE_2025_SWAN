//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#include "ChrisBrain.h"

std::shared_ptr<ParameterLink<bool>> SWNBrain::recordActivityPL = Parameters::register_parameter(
	"BRAIN_CHRIS-recordActivity", false,
	"if true, activity will be recorded to FireReport.txt");
bool SWNBrain::recordActivity;


std::shared_ptr<ParameterLink<int>> SWNBrain::initalNeuronCountPL = Parameters::register_parameter(
	"BRAIN_CHRIS-initalNeuronCount", 20,
	"number of neurons to initalize the brain with");

std::shared_ptr<ParameterLink<int>> SWNBrain::outputRulePL = Parameters::register_parameter(
	"BRAIN_CHRIS-outputRule", 1,
	"if internalUpates > 1, how is output determined?\n"
	"1: average of output over all steps\n"
	"2: final update value");

std::shared_ptr<ParameterLink<int>> SWNBrain::internalUpdatesPL = Parameters::register_parameter(
	"BRAIN_CHRIS-internalUpates", 1,
	"number of times each neuron fires per update");

// connect if ((1-distance) ^ distanceExponent) * cMatch > wireConnectionThreshold
std::shared_ptr<ParameterLink<double>> SWNBrain::wireConnectionThresholdPL =
Parameters::register_parameter(
	"BRAIN_CHRIS-wireConnectionThreshold", .5,
	"connect if ((1-distance) ^ distanceExponent) * (cMatch ^ tagExponent) > wireConnectionThreshold");
std::shared_ptr<ParameterLink<double>> SWNBrain::distanceExponentPL =
Parameters::register_parameter(
	"BRAIN_CHRIS-distanceExponent", 2.0,
	"connect if ((1-distance) ^ distanceExponent) * (cMatch ^ tagExponent) > wireConnectionThreshold\n"
	"  if -1, then distance will not effect connection");
std::shared_ptr<ParameterLink<double>> SWNBrain::cTagExponentPL =
Parameters::register_parameter(
	"BRAIN_CHRIS-cTagExponent", 1.0,
	"connect if ((1-distance) ^ distanceExponent) * (cMatch ^ tagExponent) > wireConnectionThreshold");




std::shared_ptr<ParameterLink<int>> SWNBrain::maxNeuronsPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-maxNeurons", 50,
	"mutations will not result in more then this many neurons");
std::shared_ptr<ParameterLink<int>> SWNBrain::minNeuronsPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-minNeurons", 3,
	"mutations will not result in fewer then this many neurons");

std::shared_ptr<ParameterLink<double>> SWNBrain::mutateDeleteNeuronPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateDeleteNeuron", .005,
	"per brain chance to delete a neuron");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateNewNeuronPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateNewNeuronPL", .0025,
	"per brain chance to add a random new neuron");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateCopyNeuronPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateCopyNeuron", .0025,
	"per brain chance to add a copy of an existing neuron");





std::shared_ptr<ParameterLink<double>> SWNBrain::mutateNeuronLocationPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateNeuronLocation", .001,
	"per brain chance to mutate a neurons location");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateNeuronCTagOncePL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateNeuronCTagOnce", .0025,
	"per brain chance to mutate a single bit on a neurons cTag");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateNeuronCTagRandomPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateNeuronCTagRandom", .0001,
	"per brain chance to mutate a neurons cTag to a random value");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateNeuronWTagOncePL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateNeuronWTagOnce", .0025,
	"per brain chance to mutate a single bit on a neurons wTag");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateNeuronWTagRandomPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateNeuronWTagRandom", .0001,
	"per brain chance to mutate a neurons wTag to a random value");

std::shared_ptr<ParameterLink<double>> SWNBrain::mutateOutputLocationPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateOutputLocation", .05,
	"per brain chance to mutate an output location");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateDistancePL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateDistance", .1,
	"when a mutation causes a location change, what is the max range (-1 * mutateDistance -> mutateDistance)");

// per neuron mutation rates
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateBranchCTagOncePL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateBranchCTagOnce", .001,
	"per neuron chance to flip a single bit on the cTag of one of that neurons branches");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateBranchCTagRandomPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateBranchCTagRandom", .0001,
	"per neuron chance to mutate the cTag of one of that neurons branches to a random value");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateBranchWTagOncePL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateBranchWTagOnce", .001,
	"per neuron chance to flip a single bit on the wTag of one of that neurons branches");
std::shared_ptr<ParameterLink<double>> SWNBrain::mutateBranchWTagRandomPL =
Parameters::register_parameter("BRAIN_CHRIS_MUTATIONS-mutateBranchWTagRandom", .0001,
	"per neuron chance to mutate the wTag of one of that neurons branches to a random value");



// spikey neuron parameters
std::shared_ptr<ParameterLink<double>> SWNBrain::Spikey_Neuron_Logic_Device::mutateInitalChargePL =
Parameters::register_parameter("BRAIN_CHRIS_SPIKEY_MUTATIONS-mutateInitalCharge", .005,
	"per neuron logic gate chance to mutate the Inital Charge");

std::shared_ptr<ParameterLink<double>> SWNBrain::Spikey_Neuron_Logic_Device::mutateMaxChargePL =
Parameters::register_parameter("BRAIN_CHRIS_SPIKEY_MUTATIONS-mutateMaxCharge", .005,
	"per neuron logic gate chance to mutate the max Charge");

std::shared_ptr<ParameterLink<double>> SWNBrain::Spikey_Neuron_Logic_Device::mutateMinChargePL =
Parameters::register_parameter("BRAIN_CHRIS_SPIKEY_MUTATIONS-mutateMinCharge", .005,
	"per neuron logic gate chance to mutate the min Charge");

std::shared_ptr<ParameterLink<double>> SWNBrain::Spikey_Neuron_Logic_Device::mutateThresholdPL =
Parameters::register_parameter("BRAIN_CHRIS_SPIKEY_MUTATIONS-mutateThreshold", .005,
	"per neuron logic gate chance to mutate the Threshold");

std::shared_ptr<ParameterLink<double>> SWNBrain::Spikey_Neuron_Logic_Device::mutateDeliveryChargePL =
Parameters::register_parameter("BRAIN_CHRIS_SPIKEY_MUTATIONS-mutateDeliveryCharge", .005,
	"per neuron logic gate chance to mutate the Delivery Charge");

std::shared_ptr<ParameterLink<double>> SWNBrain::Spikey_Neuron_Logic_Device::mutateDecayRatePL =
Parameters::register_parameter("BRAIN_CHRIS_SPIKEY_MUTATIONS-mutateDecayRate", .005,
	"per neuron logic gate chance to mutate the Decay Rate");

std::shared_ptr<ParameterLink<double>> SWNBrain::Spikey_Neuron_Logic_Device::mutateDefaultFirePL =
Parameters::register_parameter("BRAIN_CHRIS_SPIKEY_MUTATIONS-mutateDefaultFire", .0001,
	"per neuron logic gate chance to flip the activation behavior");



std::shared_ptr<ParameterLink<double>> SWNBrain::TuningCurve_Branch_Logic_Device::mutateTuningPointPL =
Parameters::register_parameter("BRAIN_CHRIS_BRANCH_TUNING_MUTATIONS-mutateTuningPoint", .1,
	"per branch rate to mutate one tuning curve point");



int SWNBrain::lastRecordTime = 1;

SWNBrain::SWNBrain(int _nrInNodes, int _nrOutNodes, std::shared_ptr<ParametersTable> PT_)
	: AbstractBrain(_nrInNodes, _nrOutNodes, PT_) {

	recordActivity = recordActivityPL->get(PT);
	
	initalNeuronCount = initalNeuronCountPL->get(PT);
	wireConnectionThreshold = wireConnectionThresholdPL->get(PT);
	distanceExponent = distanceExponentPL->get(PT);
	cTagExponent = cTagExponentPL->get(PT);

	outputRule = outputRulePL->get(PT);
	internalUpdates = internalUpdatesPL->get(PT);

	maxNeurons = maxNeuronsPL->get(PT);
	minNeurons = minNeuronsPL->get(PT);
	mutateDeleteNeuron = mutateDeleteNeuronPL->get(PT);
	mutateNewNeuron = mutateNewNeuronPL->get(PT);
	mutateCopyNeuron = mutateCopyNeuronPL->get(PT);

	mutateNeuronLocation = mutateNeuronLocationPL->get(PT);
	mutateNeuronCTagOnce = mutateNeuronCTagOncePL->get(PT);
	mutateNeuronCTagRandom = mutateNeuronCTagRandomPL->get(PT);
	mutateNeuronWTagOnce = mutateNeuronWTagOncePL->get(PT);
	mutateNeuronWTagRandom = mutateNeuronWTagRandomPL->get(PT);

	mutateOutputLocation = mutateOutputLocationPL->get(PT);
	mutateDistance = mutateDistancePL->get(PT);

	mutateBranchCTagOnce = mutateBranchCTagOncePL->get(PT);
	mutateBranchCTagRandom = mutateBranchCTagRandomPL->get(PT);
	mutateBranchWTagOnce = mutateBranchWTagOncePL->get(PT);
	mutateBranchWTagRandom = mutateBranchWTagRandomPL->get(PT);

	// columns to be added to ave file
	popFileColumns.clear();
	popFileColumns.push_back("neurons");
	popFileColumns.push_back("branches");
	popFileColumns.push_back("subBraches");
}

std::shared_ptr<AbstractBrain> SWNBrain::makeBrain(
	std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes) {
	std::shared_ptr<SWNBrain> newBrain = std::make_shared<SWNBrain>(nrInputValues, nrOutputValues);


	// make brain from scratch
	// make inputs and hidden neurons
	for (int i = 0; i < nrInputValues; i++) {
		newBrain->addRandomInputNeuron();
	}
	// make some neurons with branches
	for (int i = 0; i < initalNeuronCount; i++) {
		newBrain->addRandomLogicNeuron();
	}

	// make some outputs
	newBrain->addRandomOutputs();


	// 0000000000 = 0
	// 0101010101 = 341
	// 1010101010 = 682
	// 1111111111 = 1023
	// 1111100000 = 992
	// 1111000000 = 960
	// 1110000000 = 896
	// 1100000000 = 768
	// 1000000000 = 512

	// 0010010010 = 146
	// 0010010011 = 147

	// 0100100100 = 292
	// 1001001000 = 584

	/*newBrain->neurons[0].x = .5;
	newBrain->neurons[0].y = .1;
	newBrain->neurons[0].cTag = 0;
	newBrain->neurons[0].wTag = 0;


	newBrain->neurons[1].x = .5;
	newBrain->neurons[1].y = .3;
	newBrain->neurons[1].cTag = 1023;
	newBrain->neurons[1].wTag = 1023;

	newBrain->neurons[1].branches[0].cTag = 341;
	newBrain->neurons[1].branches[0].wTag = 341;
	newBrain->neurons[1].branches[1].cTag = 1023;
	newBrain->neurons[1].branches[1].wTag = 1023;
	newBrain->neurons[1].branches[2].cTag = 768;
	newBrain->neurons[1].branches[2].wTag = 768;



	newBrain->neurons[2].x = .6;
	newBrain->neurons[2].y = .3;
	newBrain->neurons[2].cTag = 1023;
	newBrain->neurons[2].wTag = 1023;

	newBrain->neurons[2].branches[0].cTag = 341;
	newBrain->neurons[2].branches[0].wTag = 1023;
	newBrain->neurons[2].branches[1].cTag = 1023;
	newBrain->neurons[2].branches[1].wTag = 341;
	newBrain->neurons[2].branches[2].cTag = 768;
	newBrain->neurons[2].branches[2].wTag = 768;
	*/

	// connect neurons and outputs
	newBrain->wireNeurons();
	newBrain->wireOutputs();

	// show brain
	/*
	for (int n = 0; n < newBrain->neurons.size(); n++) {
		std::cout << "neuron  " << n << "  ";
		if (newBrain->neurons[n].isInputNeuron) {
			std::cout << "*input*  ";
		}
		std::cout << newBrain->neurons[n].x << "," << newBrain->neurons[n].y << "  (c:" << newBrain->neurons[n].cTag << ",s:" << newBrain->neurons[n].wTag << ")" << std::endl;
		for (auto & b : newBrain->neurons[n].branches) {
			std::cout <<     "  branch (cTag:" << b.cTag << ",wTag:" << b.wTag << ")" << std::endl;
			for (int q = 0; q < b.inputConnectionIDs.size(); q++) {
				std::cout << "    " << q << "  " << b.inputConnectionIDs[q] << ", " << b.inputWeights[q] << std::endl;
			}
		}
	}
	*/

	// graph brain
	
	//std::cout << "---------------------------------------------------" << std::endl;
	//std::cout << "Brain:" << std::endl;
	//std::cout << "---------------------------------------------------" << std::endl;
	//newBrain->graphBrain();
	//std::cout << "---------------------------------------------------" << std::endl;
	
	return newBrain;
}

std::shared_ptr<AbstractBrain> SWNBrain::makeCopy(std::shared_ptr<ParametersTable> PT_) {
	auto newBrain = std::make_shared<SWNBrain>(nrInputValues, nrOutputValues, PT_);
	for (auto n : neurons) {
		newBrain->neurons.push_back(n);
	}
	for (int i = 0; i < outputNeuronIDs.size(); i++) {
		newBrain->outputNeuronIDs.push_back(outputNeuronIDs[i]);
		newBrain->outputLocations.push_back(outputLocations[i]);
	}

	if (newBrain->neurons.size() > 0) {
		newBrain->wireNeurons();
		newBrain->wireOutputs();
	}

	return newBrain;
}

std::shared_ptr<AbstractBrain>
SWNBrain::makeBrainFrom(std::shared_ptr<AbstractBrain> parent, std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes) {

	auto parentBrain = std::dynamic_pointer_cast<SWNBrain>(parent);
	auto newBrain = std::dynamic_pointer_cast<SWNBrain>(parent->makeCopy(PT));

	bool dbg = false; // debug
	bool changed = false; // did this brain change?
	int howManyMutations;
	
	if (dbg) std::cout << std::endl;

	// do we delete neurons?
	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateDeleteNeuron);
	for (int i = 0; i < howManyMutations; i++) {
		if (newBrain->neurons.size() > ((size_t)minNeurons + (size_t)newBrain->nrInputValues)) {
			if (dbg) std::cout << "DEL  ";
			changed = true;
			int which = Random::getIndex(newBrain->neurons.size() - newBrain->nrInputValues) + newBrain->nrInputValues;
			newBrain->neurons[which] = newBrain->neurons.back();
			newBrain->neurons.pop_back();
		}
	}




	// do we copy any existing neurons?
	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateCopyNeuron);
	for (int i = 0; i < howManyMutations; i++) {
		if (newBrain->neurons.size() < ((size_t)maxNeurons + (size_t)newBrain->nrInputValues)) {
			if (dbg) std::cout << "COPY  ";
			changed = true;
			newBrain->neurons.push_back(newBrain->neurons[Random::getIndex(newBrain->neurons.size()-(size_t)nrInputValues)+ (size_t)nrInputValues]);
		}
	}





	// do we create new random neurons?
	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateNewNeuron);
	for (int i = 0; i < howManyMutations; i++) {
		if (newBrain->neurons.size() < ((size_t)maxNeurons + (size_t)newBrain->nrInputValues)) {
			if (dbg) std::cout << "NEW  ";
			changed = true;
			newBrain->addRandomLogicNeuron();
		}
	}

	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateNeuronLocation);
	for (int i = 0; i < howManyMutations; i++) {
		if (dbg) std::cout << "N-LOC ";
		changed = true;
		auto &n = newBrain->neurons[Random::getIndex(newBrain->neurons.size())];
		n.x = std::max(std::min(n.x + std::pow(Random::getDouble(-1 * mutateDistance, mutateDistance), 1), 1.0), 0.0);
		n.y = std::max(std::min(n.y + std::pow(Random::getDouble(-1 * mutateDistance, mutateDistance), 1), 1.0), 0.0);
	}
	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateNeuronCTagOnce);
	for (int i = 0; i < howManyMutations; i++) {
		if (dbg) std::cout << "N-CTAG-1  ";
		changed = true;
		auto &n = newBrain->neurons[Random::getIndex(newBrain->neurons.size())];
		n.cTag.flip(Random::getIndex(tagSize));
	}
	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateNeuronCTagRandom);
	for (int i = 0; i < howManyMutations; i++) {
		if (dbg) std::cout << "N-CTAG-R  ";
		changed = true;
		auto &n = newBrain->neurons[Random::getIndex(newBrain->neurons.size())];
		n.cTag = Random::getIndex(std::pow(2, tagSize));
	}
	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateNeuronWTagOnce);
	for (int i = 0; i < howManyMutations; i++) {
		if (dbg) std::cout << "N-WTAG-1  ";
		changed = true;
		auto &n = newBrain->neurons[Random::getIndex(newBrain->neurons.size())];
		n.wTag.flip(Random::getIndex(tagSize));
	}
	howManyMutations = Random::getBinomial(newBrain->neurons.size(), mutateNeuronWTagRandom);
	for (int i = 0; i < howManyMutations; i++) {
		if (dbg) std::cout << "N-WTAG-R  ";
		changed = true;
		auto &n = newBrain->neurons[Random::getIndex(newBrain->neurons.size())];
		n.wTag = Random::getIndex(std::pow(2, tagSize));
	}

	for (auto& n : newBrain->neurons) {
		if (!n.isInputNeuron){
			howManyMutations = Random::getBinomial(n.branches.size(), mutateBranchCTagOnce);
			for (int i = 0; i < howManyMutations; i++) {
				if (dbg) std::cout << "B-CTAG-1  ";
				changed = true;
				auto& b = n.branches[Random::getIndex(n.branches.size())];
				b.cTag.flip(Random::getIndex(tagSize));
			}
			howManyMutations = Random::getBinomial(n.branches.size(), mutateBranchCTagRandom);
			for (int i = 0; i < howManyMutations; i++) {
				if (dbg) std::cout << "B-CTAG-R  ";
				changed = true;
				auto& b = n.branches[Random::getIndex(n.branches.size())];
				b.cTag = Random::getIndex(std::pow(2, tagSize));
			}
			howManyMutations = Random::getBinomial(n.branches.size(), mutateBranchWTagOnce);
			for (int i = 0; i < howManyMutations; i++) {
				if (dbg) std::cout << "B-WTAG-1  ";
				changed = true;
				auto& b = n.branches[Random::getIndex(n.branches.size())];
				b.wTag.flip(Random::getIndex(tagSize));
			}
			howManyMutations = Random::getBinomial(n.branches.size(), mutateBranchWTagRandom);
			for (int i = 0; i < howManyMutations; i++) {
				if (dbg) std::cout << "B-WTAG-R  ";
				changed = true;
				auto& b = n.branches[Random::getIndex(n.branches.size())];
				b.wTag = Random::getIndex(std::pow(2, tagSize));
			}
			n.logic->mutate();
		}
		for (auto & b : n.branches) {
			b.logic->mutate(); // call mutate on each branch logic
		}
	}
	howManyMutations = Random::getBinomial(newBrain->nrOutputValues, mutateOutputLocation);
	for (int i = 0; i < howManyMutations; i++) {
		if (dbg) std::cout << "OUT-LOC  ";
		changed = true;
		int whichOutput = Random::getIndex(newBrain->nrOutputValues);
		outputLocations[whichOutput].first = std::max(std::min(outputLocations[whichOutput].first + std::pow(Random::getDouble(-1 * mutateDistance, mutateDistance), 1), 1.0), 0.0);
		outputLocations[whichOutput].second = std::max(std::min(outputLocations[whichOutput].second + std::pow(Random::getDouble(-1 * mutateDistance, mutateDistance), 1), 1.0), 0.0);
	}
			
	newBrain->wireNeurons();
	newBrain->wireOutputs();

	return newBrain;
}

std::shared_ptr<AbstractBrain> SWNBrain::makeBrainFromMany(std::vector<std::shared_ptr<AbstractBrain>> parents,
	std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes) {
	auto parentBrain = std::dynamic_pointer_cast<SWNBrain>(parents[0]);
	auto newBrain = std::dynamic_pointer_cast<SWNBrain>(parents[0]->makeCopy(PT));

	return newBrain;
}


void SWNBrain::resetBrain()
{
	AbstractBrain::resetBrain();
	for (int i = 0; i < nrOutputValues; i++) {
		outputValues[i] = 0;
	}
	for (auto& n : neurons) {
		n.reset();
	}
}


void SWNBrain::update() {
	if (recordActivity) {
		//FileManager::writeToFile("FireReport.csv", std::to_string(nrInputValues) + "," + std::to_string(neurons.size() - nrInputValues) + ", x", "charge, threshold, status");
		FileManager::writeToFile("FireReport.txt", std::to_string(nrInputValues) + "," + std::to_string(neurons.size() - nrInputValues) + "," + std::to_string(nrOutputValues) + "," + std::to_string(internalUpdates) + ",x",
			"line endings indicate data type:\n"
			"  x : start of new brain update\n"
			"  n : neuron (active/reverse) not firing\n"
			"  f : active type neuron firing\n"
			"  r : reverse type neuron firing\n"
			"  i : input 'neuron' (just input value)\n"
			"  o : output value\n"
			"start (x) line format:\n"
			"  [# inputs],[# non-input neurons],[# outputs],[internalUpdates]\n"
			"neuron (n,f,r) line format:\n"
			"  [charge],[threshold],[active/reverse](1 indicates active),[state](n,f,r)\n"
			"input/output (i,o) line format:\n"
			"  [input value/output value]");
	}
	// record a brain from time to time
	if (Global::update%100 == 0 && lastRecordTime < Global::update ) {
		graphBrain();
		showBrain();
	}
	lastRecordTime = Global::update;

	// the first nrInputValues neurons are input neurons, load input values into these neurons
	for (int i = 0; i < nrInputValues; i++) {
		neurons[i].outputValue = inputValues[i];// (inputValues[i] != 0) ? inputValues[i] : -1;
		//if (i == nrInputValues-1){
		//	neurons[i].outputValue = 1.0;
		//}
		if (recordActivity) {
			FileManager::writeToFile("FireReport.txt", std::to_string(neurons[i].outputValue)+",i");
		}

	}

	bool db = false; // debug

	// make containers for dendrite and neurons values
	std::vector<double> dendriteValues(100, 0.0);
	std::vector<double> branchValues(100, 0.0);

	if (outputRule == 1) { // average rule
		for (int i = 0; i < nrOutputValues; i++) {
			outputValues[i] = 0;
		}
	}

	if (outputRule == 3) { // running average rule
		//for (int i = 0; i < nrOutputValues; i++) {
		//	outputValues[i] *= -1;
		//}
	}

	for (int repeat = 0; repeat < internalUpdates; repeat++) { // for number of brain updates per world update
		for (int i = nrInputValues; i < neurons.size(); i++) { // for each non-input neuron
			auto & n = neurons[i];
			for (int j = 0; j < n.branches.size(); j++) { // for each branch in the current neuron

				auto& b = n.branches[j];

				for (int k = 0; k < b.inputConnectionIDs.size(); k++) { // for each dendrite in this branch
					// dendrite = 10 * neuron value this dendrite connects to * weight for this dendrite
					dendriteValues[k] = 1.0
						* neurons[b.inputConnectionIDs[k]].outputValue
						* b.inputWeights[k];
					if (db) {
						std::cout << "neuron: " << i - nrInputValues << "  branch: " << j << "  connection: " << k <<
							" from: " << b.inputConnectionIDs[k] << " :: " <<
							"(" << neurons[b.inputConnectionIDs[k]].outputValue << " * " << b.inputWeights[k] << ") = " <<
							dendriteValues[k] << " " << std::endl;
					}
				}
				// now that we know the dendrite values, update the branchValues using the branch logic 
				
				

				branchValues[j] = n.branchWeights[j] * b.logic->update(dendriteValues, b.inputConnectionIDs.size());
				//branchValues[j] = b.logic->update(dendriteValues, b.inputConnectionIDs.size());
				
				
				
				if (db) {
					std::cout << "     post branch logic = " << branchValues[j] << "    \n";
				}
			}
			// now that we have all the branch values, run the neurons logic
			//std::cout << "   r:"<< repeat << "  neuron#:" << i << std::endl;
			n.outputValue_new = n.logic->update(branchValues, n.branches.size());
			if (0) {
				for (int j = 0; j < n.branches.size(); j++) {
					// for each branch, if the branch is active and the neuron logic is active, increse the branch weight
					if (std::abs(branchValues[j]) > .001 && std::abs(n.outputValue_new) > .001) {
						n.branchWeights[j] *= 1.01;
					}
					// if branch is active and neuron logic is not, then decrese branch weight
					else if (std::abs(branchValues[j]) > .001 && std::abs(n.outputValue_new) < .001) {
						n.branchWeights[j] *= .99;
					}
				}
			}
			if (db) {
				std::cout << " gate out: " << n.outputValue_new << std::endl;
			}
		}
		// once all the neurons have been run, update each neurons outputValues
		for (int i = nrInputValues; i < neurons.size(); i++) {
			neurons[i].outputValue = neurons[i].outputValue_new;
		}
		if (outputRule == 1 || outputRule == 3) { // average rule
			for (int i = 0; i < nrOutputValues; i++) {
				// if rule is 1 or 3, update brains output values based on current neuron states divided by number of brain updates per world update
				outputValues[i] += neurons[outputNeuronIDs[i]].outputValue/static_cast<double>(internalUpdates);
				//if (recordActivity) {
				//	FileManager::writeToFile("FireReport.txt", std::to_string(outputValues[i]) + ",o");
				//}
			}
		}
	}

	if (outputRule == 2) { // final rule
		// if rule is 2, update brain output with current neuron states
		for (int i = 0; i < nrOutputValues; i++) {
			outputValues[i] = neurons[outputNeuronIDs[i]].outputValue;
		}
	}

	if (recordActivity) {
		for (int i = 0; i < nrOutputValues; i++) {
			FileManager::writeToFile("FireReport.txt", std::to_string(outputValues[i]) + ",o");
		}
	}

	if (db) {
		std::cout << "   outs: ";
		for (int i = 0; i < nrOutputValues; i++) {
			std::cout << outputValues[i] << "  ";
		}
		std::cout << std::endl;
	}

}

std::string SWNBrain::description() { return "Chris Brain\n"; }

DataMap SWNBrain::getStats(std::string& prefix) { 
	DataMap tempDM;
	int totalBranches = 0;
	int totalSubBraches = 0;
	for (auto & n : neurons) {
		totalBranches += n.branches.size();
		for (auto & b : n.branches) {
			totalSubBraches += b.inputConnectionIDs.size();
		}
		tempDM.append(prefix+"neurons", (int)(neurons.size() - nrInputValues));
		tempDM.append(prefix+"branches", totalBranches);
		tempDM.append(prefix+"subBraches", totalSubBraches);
	}
	return tempDM;
}

void SWNBrain::initializeGenomes(
	std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes) {
	// do nothing;
}




// convert a brain into data map with data that can be saved to file so this brain can be reconstructed
// 'name' here contains the prefix that must be so that deserialize can identify relavent data
DataMap SWNBrain::serialize(std::string& name) {
	
	//showBrain();

	DataMap dataMap;
	std::stringstream ss;

	for (auto ol : outputLocations) {
		ss << ol.first << "^" << ol.second << "%";
	}
	
	for (auto neuron : neurons) {
		/*
		bool isInputNeuron; // is this an input neuron? if yes, then outvalue will be set on update and neuron is not itself run
		std::vector<DendriticArbor> branches;
		std::vector<double> branchWeights; // value mutiplier for each branch
		double x, y;                       // location of neuron (for connection tag)
		std::bitset<tagSize> cTag;         // conection tag
		std::bitset<tagSize> wTag;         // weight tag
		double outputValue = 0.0;          // value generated by neuron this update
		double outputValue_new = 0.0;      // current value in this neuron when other neurons access it

		std::shared_ptr<Abstract_Neuron_Logic_Device> logic;   // logic in neuron (this will take branch values as input and output to outputValue_new
		*/


		ss << "N" << (int)neuron.isInputNeuron << "<" << neuron.x << "<" << neuron.y << "<" <<
			neuron.cTag << "<" << neuron.wTag << "<";


		if (!neuron.isInputNeuron) {
			ss << neuron.logic->serialize().str();
			for (int i = 0; i < neuron.branches.size(); i++) {
				ss << "*" << neuron.branchWeights[i] << "|" << neuron.branches[i].cTag <<
					"|" << neuron.branches[i].wTag << "|" << neuron.branches[i].logic->serialize().str();
			}
		}
		//ss.seekp(-1, std::ios_base::end); // set write head back one char
	}
	//dataMap.set(name + "brainData", ss.str().substr(0, ss.str().size() - 1));
	dataMap.set(name + "brainData", ss.str());
	//FileManager::writeToFile("SWNB_data.csv", ss.str());
	return dataMap;
}

// given an unordered_map<string, string> of org data and PT, load data into this brain
// 'name' here contains the prefix that was used when data was being saved
void SWNBrain::deserialize(std::shared_ptr<ParametersTable> PT, std::unordered_map<std::string, std::string>& orgData, std::string& name) {
		
	neurons.clear();

	std::string data = orgData[name + "brainData"];
	std::vector<std::string> BrainStr; // a list with all the brains parts
	std::vector<std::string> allNeuronsStr; // a list with all the neurons (after input locations are removed)
	std::vector<std::string> neuronLogicStr; // a list with all the data for a single neuron
	std::vector<std::string> allBranchDataStr; // a list with all branches of a single neuron
	std::vector<std::string> branchDataStr; // a list with the data for a single branch

	std::stringstream converter; // used to convert strings to other types

	convertCSVListToVector(data, BrainStr, '%'); // first division of data, isolates each outputLocation and everything else in last string
	
	// load output locations
	outputLocations.clear();
	for (int i = 0; i < BrainStr.size() - 1; i++) {
		std::vector<double> loc;
		convertCSVListToVector(BrainStr[i], loc,'^');
		outputLocations.push_back({ loc[0], loc[1] });
	}

	convertCSVListToVector(BrainStr.back(), allNeuronsStr, 'N'); // divide remaining data into neurons
	
	// iterate over neurons in data
	// since neurons start with "N" there is a blank vector at the head, (so start at [1])
	for (int n = 1; n < allNeuronsStr.size(); n++) {
		//std::cout << " n = " << n << std::endl;
		std::vector<std::string> neuronStr; // a list with the parts a neuron
		convertCSVListToVector(allNeuronsStr[n], neuronStr, '<'); // isolate parts of neuron (inputNeuron, loc, c&w tags, logic and branches)

		if (neuronStr[0] == "1") {
			//std::cout << "   1" << std::endl;
			neurons.push_back(Neuron(true));
		}
		else {
			//std::cout << "   0" << std::endl;
			neurons.push_back(Neuron(false));
		}
		neurons.back().branches.clear();
		neurons.back().branchWeights.clear();
		converter.clear();
		converter.str(neuronStr[1]);
		converter >> neurons.back().x;
		converter.clear(); 
		converter.str(neuronStr[2]);
		converter >> neurons.back().y;
		converter.clear();
		converter.str(neuronStr[3]);
		converter >> neurons.back().cTag;
		converter.clear();
		converter.str(neuronStr[4]);
		converter >> neurons.back().wTag;
		if (neuronStr[0] == "0") { //if this is NOT an input neuron
			convertCSVListToVector(neuronStr[5], neuronLogicStr, '/');
			convertCSVListToVector(neuronLogicStr.back(), allBranchDataStr,'*');
			// load the logic
			if (neuronLogicStr[0] == "snld") { // this is a spikey neuron logic
				auto newLogic = std::make_shared<Spikey_Neuron_Logic_Device>(allBranchDataStr.size() - 1);
				newLogic->deserialize(neuronLogicStr);
				neurons.back().logic = newLogic;
			}
			else {
				// other neuron logic loading here
			}

			// now deal with the branches
			int branchCount = 0;
			for (int b = 1; b < allBranchDataStr.size(); b++) {
				neurons.back().branches.push_back(DendriticArbor());
				neurons.back().branches.back().inputConnectionIDs.clear();
				neurons.back().branches.back().inputWeights.clear();

				convertCSVListToVector(allBranchDataStr[b], branchDataStr, '|');
				converter.clear();
				converter.str(branchDataStr[0]);
				double w;
				converter >> w;
				neurons.back().branchWeights.push_back(w);
				converter.clear();
				converter.str(branchDataStr[1]);
				converter >> neurons.back().branches.back().cTag;
				converter.clear();
				converter.str(branchDataStr[2]);
				converter >> neurons.back().branches.back().wTag;


				if (branchDataStr[3] == "ACC_branchLogic") {
					neurons.back().branches.back().logic = std::make_shared<Accumulator_Branch_Logic_Device>(Accumulator_Branch_Logic_Device());
				}
				if (branchDataStr[3] == "TC_branchLogic") {
					neurons.back().branches.back().logic = std::make_shared<TuningCurve_Branch_Logic_Device>(TuningCurve_Branch_Logic_Device());
					neurons.back().branches.back().logic->deserialize(branchDataStr);
				}
			}
		}
	}
	wireNeurons();
	wireOutputs();	
	
	resetBrain();
	showBrain();

	graphBrain();

	auto test_serialize = serialize(name);
	test_serialize.openAndWriteToFile("SWNB_data.csv");
	std::cout << "saved a copy to SWNB_data.csv" << std::endl;
}

