//  MABE is a product of The Hintze Lab @ MSU
//     for general research information:
//         hintzelab.msu.edu
//     for MABE documentation:
//         github.com/Hintzelab/MABE/wiki
//
//  Copyright (c) 2015 Michigan State University. All rights reserved.
//     to view the full license, visit:
//         github.com/Hintzelab/MABE/wiki/License

#pragma once

#include <cmath>
#include <memory>
#include <iostream>
#include <set>
#include <vector>
#include <bitset>

#include "../../Genome/AbstractGenome.h"

#include "../../Utilities/Random.h"

#include "../AbstractBrain.h"

class SWNBrain : public AbstractBrain {
public:

	static int lastRecordTime;

	//static std::shared_ptr<ParameterLink<bool>> useActionMapPL;
	///////////////////////////////////////////////////////////////////////////////////
	// Parameters /////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	static const int tagSize = 8;

	static std::shared_ptr<ParameterLink<bool>> recordActivityPL;
	static bool recordActivity;

	static std::shared_ptr<ParameterLink<int>> outputRulePL;
	int outputRule;
	static std::shared_ptr<ParameterLink<int>> internalUpdatesPL;
	int internalUpdates;

	static std::shared_ptr<ParameterLink<int>> initalNeuronCountPL;
	int initalNeuronCount;

	static std::shared_ptr<ParameterLink<int>> maxNeuronsPL;
	int maxNeurons;
	static std::shared_ptr<ParameterLink<int>> minNeuronsPL;
	int minNeurons;

	// connect if ((1-distance) ^ distanceExponent) * cMatch > wireConnectionThreshold
	static std::shared_ptr<ParameterLink<double>> wireConnectionThresholdPL;
	double wireConnectionThreshold;
	static std::shared_ptr<ParameterLink<double>> distanceExponentPL;
	double distanceExponent;
	static std::shared_ptr<ParameterLink<double>> cTagExponentPL;
	double cTagExponent;

	// per brain mutation rates
	static std::shared_ptr<ParameterLink<double>> mutateDeleteNeuronPL;
	double mutateDeleteNeuron;
	static std::shared_ptr<ParameterLink<double>> mutateNewNeuronPL;
	double mutateNewNeuron;
	static std::shared_ptr<ParameterLink<double>> mutateCopyNeuronPL;
	double mutateCopyNeuron;


	static std::shared_ptr<ParameterLink<double>> mutateNeuronLocationPL;
	double mutateNeuronLocation;
	static std::shared_ptr<ParameterLink<double>> mutateNeuronCTagOncePL;
	double mutateNeuronCTagOnce;
	static std::shared_ptr<ParameterLink<double>> mutateNeuronCTagRandomPL;
	double mutateNeuronCTagRandom;
	static std::shared_ptr<ParameterLink<double>> mutateNeuronWTagOncePL;
	double mutateNeuronWTagOnce;
	static std::shared_ptr<ParameterLink<double>> mutateNeuronWTagRandomPL;
	double mutateNeuronWTagRandom;

	static std::shared_ptr<ParameterLink<double>> mutateOutputLocationPL;
	double mutateOutputLocation;
	static std::shared_ptr<ParameterLink<double>> mutateDistancePL;
	double mutateDistance;

	// per neuron mutation rates
	static std::shared_ptr<ParameterLink<double>> mutateBranchCTagOncePL;
	double mutateBranchCTagOnce;
	static std::shared_ptr<ParameterLink<double>> mutateBranchCTagRandomPL;
	double mutateBranchCTagRandom;
	static std::shared_ptr<ParameterLink<double>> mutateBranchWTagOncePL;
	double mutateBranchWTagOnce;
	static std::shared_ptr<ParameterLink<double>> mutateBranchWTagRandomPL;
	double mutateBranchWTagRandom;


	///////////////////////////////////////////////////////////////////////////////////
	// Utilites ///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	// return # of matching bits / size of bitset
	template <class T>
	double bitMatch(T a, T b) {
		return (static_cast<double>(((a ^ b)).count()) / static_cast<double>(a.size()));
	}

	// return distance between x1, y1 and x2, y2
	double dist(double x1, double x2, double y1, double y2) {
		return (std::sqrt(std::pow(x1 - x2, 2.0) + std::pow(y1 - y2, 2.0)));
	}

	///////////////////////////////////////////////////////////////////////////////////
	// Dedrite_Logic_Devices //////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	// Abstract_Dedrite_Logic_Device must be able to work with any number of inputs
	// because we don't know how many inputs a DendriticArbor may have.
	class Abstract_Branch_Logic_Device {
	public:
		virtual double update(const std::vector<double> & inputValues,int nrValues) {
			std::cout << "attempt to use a Dedrite_Logic_Device without a defined update function.\n  exiting..." << std::endl;
			exit(1);
		}
		virtual std::stringstream serialize() = 0;
		virtual void deserialize(std::vector<std::string> branchData) = 0;

		virtual std::shared_ptr< Abstract_Branch_Logic_Device> makeCopy() = 0;

		virtual void mutate() {};

	};

	// return ave of inputValues
	class Accumulator_Branch_Logic_Device : public Abstract_Branch_Logic_Device {
	public:
		double update(const std::vector<double>& inputValues, int nrValues) override {
			if (nrValues == 0 || inputValues.size() == 0) {
				return 0.0;
			}
			//return(std::accumulate(inputValues.begin(), inputValues.begin() + nrValues, 0.0) / static_cast<double>(nrValues));
			return(std::accumulate(inputValues.begin(), inputValues.begin() + nrValues, 0.0));
		}
		std::stringstream serialize() override {
			std::stringstream ss;
			ss << "ACC_branchLogic";
			return ss;
		}
		void deserialize(std::vector<std::string> branchData) override {
			// do nothing, this logic type has no mutable aspects
		}

		virtual std::shared_ptr< Abstract_Branch_Logic_Device> makeCopy() {
			return(std::make_shared<Accumulator_Branch_Logic_Device>(Accumulator_Branch_Logic_Device()));
		}

		// using the default mutate since this has no mutable parts

	};

	// return ave of inputValues
	class TuningCurve_Branch_Logic_Device : public Abstract_Branch_Logic_Device {
	public:

		static std::shared_ptr<ParameterLink<double>> mutateTuningPointPL;
		std::vector<double> tuningCurve;

		TuningCurve_Branch_Logic_Device() {
			//tuningCurve = { Random::getDouble(1.0),Random::getDouble(1.0),Random::getDouble(1.0),Random::getDouble(1.0),Random::getDouble(1.0) };
			tuningCurve = { 0.0,0.25,0.5,0.75,1.0 };
		}

		TuningCurve_Branch_Logic_Device(std::vector<double> _tuningCurve) :tuningCurve(_tuningCurve){
		}



		double update(const std::vector<double>& inputValues, int nrValues) override {
			if (nrValues == 0 || inputValues.size() == 0) {
				return 0.0;
			}
			//return(std::accumulate(inputValues.begin(), inputValues.begin() + nrValues, 0.0) / static_cast<double>(nrValues));
			double accumulant = std::accumulate(inputValues.begin(), inputValues.begin() + nrValues, 0.0);
			if (accumulant <= 0) { return(tuningCurve[0]); }
			if (accumulant > 0 && accumulant <= .25) {
				double x = 4 * accumulant;
				return((1 - x) * tuningCurve[0] + x * tuningCurve[1]);
			}
			if (accumulant > .25 && accumulant <= .5) {
				double x = 4 * (accumulant - .25);
				return((1 - x) * tuningCurve[1] + x * tuningCurve[2]);
			}
			if (accumulant > .5 && accumulant <= .75) {
				double x = 4 * (accumulant - .5);
				return((1 - x) * tuningCurve[2] + x * tuningCurve[3]);
			}
			if (accumulant > .75 && accumulant <= 1) {
				double x = 4 * (accumulant - .75);
				return((1 - x) * tuningCurve[3] + x * tuningCurve[4]);
			}
			else { return(tuningCurve[4]); } // (accumulant > 1)
		}

		std::stringstream serialize() override {
			std::stringstream ss;
			ss << "TC_branchLogic|" << tuningCurve[0] << "|" << tuningCurve[1] << "|" << 
				tuningCurve[2] << "|" << tuningCurve[3] << "|" << tuningCurve[4];
			return ss;
		}

		void deserialize(std::vector<std::string> branchData) override {
			std::stringstream converter; // used to convert strings to other types
			converter.clear();
			converter.str(branchData[4]);
			converter >> tuningCurve[0];
			converter.clear();
			converter.str(branchData[5]);
			converter >> tuningCurve[1];
			converter.clear();
			converter.str(branchData[6]);
			converter >> tuningCurve[2];
			converter.clear();
			converter.str(branchData[7]);
			converter >> tuningCurve[3];
			converter.clear();
			converter.str(branchData[8]);
			converter >> tuningCurve[4];
			std::cout << "deserialize TC_branchLogic: got... ";
			std::cout << tuningCurve[0] << "," << tuningCurve[1] << "," << tuningCurve[2] << "," <<
				tuningCurve[3] << "," << tuningCurve[4] << std::endl;
		}

		virtual std::shared_ptr< Abstract_Branch_Logic_Device> makeCopy() {
			return(std::make_shared<TuningCurve_Branch_Logic_Device>(TuningCurve_Branch_Logic_Device(tuningCurve)));
		}

		virtual void mutate() {
			if (Random::P(mutateTuningPointPL->get())) {
				tuningCurve[Random::getIndex(tuningCurve.size())] += Random::getDouble(-.25, .25);
			}
		};


	};


	///////////////////////////////////////////////////////////////////////////////////
	// Neuorn_Logic_Devices ///////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	// Abstract_Neuron_Logic_Device, defines nrInputValues. nrInputValues
	// allows a) the logic in this device to assume the exact number of inputs
	// that will be provided and is used to estabish the number of branches in
	// a neuron
	class Abstract_Neuron_Logic_Device {
	public:
		int nrInputValues; // number of inputs to this I/O device
		Abstract_Neuron_Logic_Device() {
			nrInputValues = 0;
			std::cout << "attempt to construct an Abstract_Neuron_Logic.\n  exiting..." << std::endl;
			exit(1);
		};

		virtual std::shared_ptr< Abstract_Neuron_Logic_Device> makeCopy() {
			std::cout << "attempt to use a Neuron_Logic_Device makeCopy().\n  exiting..." << std::endl;
			exit(1);
		}

		Abstract_Neuron_Logic_Device(int defaultNrNeuronLogicInputs_) : nrInputValues(defaultNrNeuronLogicInputs_) {};

		virtual double update(const std::vector<double> & inputValues, int nrValues) {
			std::cout << "attempt to use a Neuron_Logic_Device without a defined update function.\n  exiting..." << std::endl;
			exit(1);
		}
		virtual void reset() {
			std::cout << "attempt to reset an Abstract_Neuron_Logic.\n  exiting..." << std::endl;
			exit(1);
		};
		virtual void mutate() {
			std::cout << "attempt to mutate an Abstract_Neuron_Logic.\n  exiting..." << std::endl;
			exit(1);
		};

		virtual std::stringstream serialize() = 0;
		virtual void deserialize(std::vector<std::string> neuronLogicStr) = 0;
	};


	class Spikey_Neuron_Logic_Device : public Abstract_Neuron_Logic_Device {
	public:

		static std::shared_ptr<ParameterLink<double>> mutateInitalChargePL;
		static std::shared_ptr<ParameterLink<double>> mutateMaxChargePL;
		static std::shared_ptr<ParameterLink<double>> mutateMinChargePL;
		static std::shared_ptr<ParameterLink<double>> mutateThresholdPL;
		static std::shared_ptr<ParameterLink<double>> mutateDeliveryChargePL;
		static std::shared_ptr<ParameterLink<double>> mutateDecayRatePL;
		static std::shared_ptr<ParameterLink<double>> mutateDefaultFirePL;

		double mutateInitalCharge;
		double mutateMaxCharge;
		double mutateMinCharge;
		double mutateThreshold;
		double mutateDeliveryCharge;
		double mutateDecayRate;
		double mutateDefaultFire;

		double currentChargeInit = 0;
		double maxChargeInit = 10;
		double minChargeInit = -10;
		double thresholdInit = .2;
		double deliveryChargeInit = 1;
		double decayRateInit = .1;

		double currentCharge = 0;
		double maxCharge = 10;
		double minCharge = -10;
		double threshold = .2;
		double deliveryCharge = 1;
		double decayRate = .1;
		bool defalutFire = true;

		std::shared_ptr< Abstract_Neuron_Logic_Device> makeCopy() {
			auto newLogic = std::make_shared<Spikey_Neuron_Logic_Device>(nrInputValues);
			newLogic->mutateInitalCharge;
			newLogic->mutateMaxCharge = mutateMaxCharge;
			newLogic->mutateMinCharge = mutateMinCharge;
			newLogic->mutateThreshold = mutateThreshold;
			newLogic->mutateDeliveryCharge = mutateDeliveryCharge;
			newLogic->mutateDecayRate = mutateDecayRate;
			newLogic->mutateDefaultFire = mutateDefaultFire;

			newLogic->currentChargeInit = currentChargeInit;
			newLogic->maxChargeInit = maxChargeInit;
			newLogic->minChargeInit = minChargeInit;
			newLogic->thresholdInit = thresholdInit;
			newLogic->deliveryChargeInit = deliveryChargeInit;
			newLogic->decayRateInit = decayRateInit;
			newLogic->defalutFire = defalutFire;
			newLogic->reset(); // set current values
			return(newLogic);
		}

		Spikey_Neuron_Logic_Device(int defaultNrNeuronLogicInputs_) : Abstract_Neuron_Logic_Device(defaultNrNeuronLogicInputs_) {
			mutateInitalCharge = mutateInitalChargePL->get();
			mutateMaxCharge = mutateMaxChargePL->get();
			mutateMinCharge = mutateMinChargePL->get();
			mutateThreshold = mutateThresholdPL->get();
			mutateDeliveryCharge = mutateDeliveryChargePL->get();
			mutateDecayRate = mutateDecayRatePL->get();
			mutateDefaultFire = mutateDefaultFirePL->get();

			currentChargeInit = Random::getDouble(-1, 0);
			maxChargeInit = 10;// Random::getDouble(1, 10);
			minChargeInit = -10;// Random::getDouble(-1, -10);
			thresholdInit = Random::getDouble(-1, 1);
			deliveryChargeInit = 1; //Random::getDouble(-2, 2);
			decayRateInit = Random::getDouble(0, 1);
			defalutFire = 1; //Random::getIndex(2);
			reset(); // set current values
		}

		double update(const std::vector<double>& inputValues, int nrValues) override  {
			
			bool db = false; // debug
			
			if (db) {
				std::cout << "           Current Charge:" << currentCharge;
			}

			//currentCharge increases by ave of inputs
			currentCharge+=(std::accumulate(inputValues.begin(), inputValues.begin() + nrValues, 0.0));
			
			if (db) {
				std::cout << " -> " << currentCharge << "  Threshold: " << threshold;
			}
			
			std::string outStr;
			if (recordActivity) {
				//std::cout << currentCharge << " : " << threshold << std::endl;
				outStr = std::to_string(currentCharge) + "," + std::to_string(threshold) + ",";
			}
			if (defalutFire) {
				if (currentCharge > threshold) {
					if (recordActivity) {
						FileManager::writeToFile("FireReport.txt", outStr + std::to_string(defalutFire) + ",f");
					}
					//currentCharge = currentChargeInit;
					currentCharge = thresholdInit + currentChargeInit;

					if (db) {
						std::cout << "   FIRE" << std::endl;
					}

					return deliveryCharge;
				}
			}
			else {
				if (currentCharge < threshold) {
					if (recordActivity) {
						FileManager::writeToFile("FireReport.txt", outStr + std::to_string(defalutFire) + ",r");
					}
					//currentCharge = currentChargeInit;
					currentCharge = thresholdInit + currentChargeInit;
					return deliveryCharge;
				}
			}
			if (recordActivity) {
				FileManager::writeToFile("FireReport.txt", outStr + std::to_string(defalutFire) + ",n");
			}
			currentCharge *= decayRate;
			currentCharge = std::min(maxCharge,std::max(minCharge, currentCharge));
			
			if (db) {
				std::cout << "   -- new current chrage (DR:" << decayRate << "): " << currentCharge << std::endl;
			}
			return 0.0;
		}

		void reset() override {
			//std::cout << "In RESET dc: " << deliveryCharge;
			currentCharge = thresholdInit + currentChargeInit;
			maxCharge = maxChargeInit;
			minCharge = minChargeInit;
			threshold = thresholdInit;
			deliveryCharge = deliveryChargeInit;
			//std::cout << deliveryCharge;
			decayRate = decayRateInit;
		}
		void mutate() override {
			bool dbg = false;// true;

			if (Random::P(mutateInitalCharge)) {
				if (dbg) std::cout << "G-CH  ";
				currentChargeInit += Random::getDouble(-.1, .1);
				currentChargeInit = std::max(-1.0, std::min(0.0, currentChargeInit));
			}
			//if (Random::P(mutateMaxCharge)) {
			//	if (dbg) std::cout << "G-MAX  ";
			//	maxChargeInit += Random::getDouble(-.1, .1);
			//}
			//if (Random::P(mutateMinCharge)) {
			//	if (dbg) std::cout << "G-MIN  ";
			//	minChargeInit += Random::getDouble(-.1, .1);
			//}
			if (Random::P(mutateThreshold)) {
				if (dbg) std::cout << "G-THRESH  ";
				thresholdInit += Random::getDouble(-.1, .1);
				thresholdInit = std::max(0.0, std::min(1.0, thresholdInit));
			}
			//if (Random::P(mutateDeliveryCharge)) {
			//	if (dbg) std::cout << "G-DCH " << deliveryChargeInit;
			//	deliveryChargeInit += Random::getDouble(-.1, .1);
			//	deliveryChargeInit = std::max(.25, std::min(2.0, deliveryChargeInit));
			//	if (dbg) std::cout << " " << deliveryChargeInit << "  ";
			//}
			if (Random::P(mutateDecayRate)) {
				if (dbg) std::cout << "G-DERATE  ";
				decayRateInit += Random::getDouble(-.1, .1);
				decayRateInit = std::max(0.0, std::min(1.0, decayRateInit));
			}
			//if (Random::P(mutateDefaultFire)) {
			//	if (dbg) std::cout << "G-FIRE  ";
			//	defalutFire += Random::getIndex(2);
			//}
		}

		std::stringstream serialize() override {
			std::stringstream ss;
			ss << "snld/" << currentChargeInit << "/" << maxChargeInit << "/" << minChargeInit << "/" <<
				thresholdInit << "/" << deliveryChargeInit << "/" << decayRateInit << "/" <<
				defalutFire << "/";
			return ss;
		}

		void deserialize(std::vector<std::string> neuronLogicStr) override {
			std::stringstream converter; // used to convert strings to other types
			converter.clear();
			converter.str(neuronLogicStr[1]);
			converter >> currentChargeInit;
			converter.clear();
			converter.str(neuronLogicStr[2]);
			converter >> maxChargeInit;
			converter.clear();
			converter.str(neuronLogicStr[3]);
			converter >> minChargeInit;
			converter.clear();
			converter.str(neuronLogicStr[4]);
			converter >> thresholdInit;
			converter.clear();
			converter.str(neuronLogicStr[5]);
			converter >> deliveryChargeInit;
			converter.clear();
			converter.str(neuronLogicStr[6]);
			converter >> decayRateInit;
			converter.clear();
			converter.str(neuronLogicStr[7]);
			converter >> defalutFire;
		}


	};


	///////////////////////////////////////////////////////////////////////////////////
	// DendriticArbor /////////////////////////////////////////////////////////////////
	// each neuron has a list of DendriticArbors
	// each DendriticArbor is a collection of dedrites (i.e. "wires") into a neuron from other neurons
	// (including input neurons) that has a cTag and wTag
	//   cTag (connection tag) determins which other neurons the dedrites in this Arbor will connect to
	//   wTag (weight tag) determins how strong each dendrite connection is (i.e. how much of the singal
	//        generated by a neuron comes though this wire.)
	//   cTax and wTag are sued to populate inputConnectionIDs and inputWeights
	// each DendriticArbor also has a logic and output. The connected neurons * dendrite tag are delivered
	// to the logic and the result it placed in outputValue
	///////////////////////////////////////////////////////////////////////////////////

	class DendriticArbor {
	public:
		std::bitset<tagSize> cTag; // conection tag
		std::bitset<tagSize> wTag; // weight tag
		double outputValue;

		// neuronConectionsIDs has list of indexes into brains neuron vector
		std::vector<int> inputConnectionIDs;
		// inputWeights contains weights for each entry in inputConnectionIDs
		std::vector<double> inputWeights;

		std::shared_ptr<Abstract_Branch_Logic_Device> logic;

		DendriticArbor(const DendriticArbor& da) {                              // copy constructor
			cTag = da.cTag;
			wTag = da.wTag;
			//outputValue = da.outputValue;
			//inputConnectionIDs = da.inputConnectionIDs;
			//inputWeights = da.inputWeights;
			logic = da.logic->makeCopy();
		}

		DendriticArbor() {
			cTag = Random::getIndex(std::pow(2, tagSize));
			wTag = Random::getIndex(std::pow(2, tagSize));
			
			logic = std::make_shared<Accumulator_Branch_Logic_Device>(Accumulator_Branch_Logic_Device());
			
			//bool disable_TuningCurve_Branch_Logic = true;

			//if (disable_TuningCurve_Branch_Logic) {
				// this will make only Accumulator_Branch_Logic
			//}
			//else {
				// this will make Accumulator_Branch_Logic and TuningCurve_Branch_Logic
				//if (Random::P(.5)) {
					//logic = std::make_shared<Accumulator_Branch_Logic_Device>(Accumulator_Branch_Logic_Device());
				//}
				//else {
					//logic = std::make_shared<TuningCurve_Branch_Logic_Device>(TuningCurve_Branch_Logic_Device());
				//}
			//}
		}

		//DendriticArbor(const DendriticArbor& b) { // old copy constructor - before branch logic had mutable attributes
		//	cTag = b.cTag;
		//	wTag = b.wTag;
		//	// this needs to be a logic copy! (also see deserialize)
		//	logic = std::make_shared<Accumulator_Branch_Logic_Device>(Accumulator_Branch_Logic_Device());
		//}
	};

	///////////////////////////////////////////////////////////////////////////////////
	// Neurons ////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	class Neuron {
	public:
		bool isInputNeuron; // is this an input neuron? if yes, then outvalue will be set on update and neuron is not itself run
		std::vector<DendriticArbor> branches;
		std::vector<double> branchWeights; // value mutiplier for each branch
		double x, y;                       // location of neuron (for connection tag)
		std::bitset<tagSize> cTag;         // conection tag
		std::bitset<tagSize> wTag;         // weight tag
		double outputValue = 0.0;          // value generated by neuron this update
		double outputValue_new = 0.0;      // current value in this neuron when other neurons access it

		std::shared_ptr<Abstract_Neuron_Logic_Device> logic;   // logic in neuron (this will take branch values as input and output to outputValue_new

		Neuron() = delete;

		Neuron(const Neuron& n) {                              // copy constructor
			isInputNeuron = n.isInputNeuron;
			x = n.x;
			y = n.y;
			cTag = n.cTag;
			wTag = n.wTag;
			if (!isInputNeuron) {
				logic = n.logic->makeCopy();
				for (auto b : n.branches) {
					branches.push_back(b);
					branchWeights.push_back(1.0);
				}
			}
		}

		// construtor to make a new neuon, takes isInputNeuron, and randomizes the rest
		// if true, isInputNeuron set true, no logic and no branches
		// if false, isInputNeuron set false, neuron will have logic and branches
		Neuron(bool isInputNeuron_) : isInputNeuron(isInputNeuron_) {
			x = Random::getDouble(0.0, 1.0);
			y = Random::getDouble(0.0, 1.0);
			cTag = Random::getIndex(std::pow(2, tagSize));
			wTag = Random::getIndex(std::pow(2, tagSize));
			if (!isInputNeuron) {
				logic = std::shared_ptr<Spikey_Neuron_Logic_Device>(new Spikey_Neuron_Logic_Device(Random::getInt(1,5))); // 1 to 6 inputs to logic (i.e. 3 DendriticArbors) 
				// make some branches
				for (int j = 0; j < logic->nrInputValues; j++) {
					branches.push_back(DendriticArbor()); // for each branche add a randomly constructed DendriticArbor
					branchWeights.push_back(1.0);         // set all inital branch weights to 1.0
				}
			}
		}

		void reset() { // to reset a neuron, reset it's output value(s) and call reset on logic (incase it has state)
			           // branches (at least now) are stateless
			outputValue = 0.0;
			outputValue_new = 0.0;
			if (!isInputNeuron) {
				logic->reset();
				// reset branch weights
				for (int j = 0; j < logic->nrInputValues; j++) {
					branchWeights[j] = 1.0;         // set all inital branch weights to 1.0
				}
			}
		}

	};

	////////////////////////////////////////////////////////////////////////////////////////////
	// Parts ///////////////////////////////////////////////////////////////////////////////////
	// a brain is made up of:
	//   1) neurons (input neurons are included)
	//   2) output locations (these will wire to the closest neuron used on construction)
	//   3) outputNeuronIDs - the ID of the neuron that each output wires too (used on update)
	////////////////////////////////////////////////////////////////////////////////////////////

	std::vector <Neuron> neurons;
	std::vector <std::pair<double, double>> outputLocations;
	std::vector <int> outputNeuronIDs;

	///////////////////////////////////////////////////////////////////////////////////
	// Neruon/Dedrite and Ouput Initialization Functions //////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	void addRandomInputNeuron() {
		auto what = Neuron(true);
		neurons.emplace_back(Neuron(true)); // call neuron constrctor to create an input neuron
	}

	void addRandomLogicNeuron() {
		neurons.emplace_back(Neuron(false)); // call neuron constrctor to create a logic neuron
	}

	void addRandomOutputs() {
		for (int i = 0; i < nrOutputValues; i++) {
			outputLocations.push_back({ Random::getDouble(0.0,1.0),Random::getDouble(0.0,1.0) });
		}
	}

	///////////////////////////////////////////////////////////////////////////////////
	// Wiring Functions ///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	// iterate over all dendrite branches and compair to all possible inputs (input
	// and logic neurons). Establish conections and conection weights
	void wireNeurons() {
		// go over all neurons (not input neurons) and create synaptic
		// connections for each neuron, skip input neurons
		for (int nID = nrInputValues; nID < neurons.size(); nID++) {
			// for each branch in each neuron
			auto& n = neurons[nID];

			for (int bID = 0; bID < n.branches.size(); bID++) {
				auto& b = n.branches[bID];
				b.inputConnectionIDs.clear();
				b.inputWeights.clear();
				// compair each branch to each neuron it's not part of
				for (int sourceID = 0; sourceID < neurons.size(); sourceID++) {
					if (sourceID != nID) {
						auto& source = neurons[sourceID];
						// see if we make a connection from b to source
						// if so, determin strength, and then add nID, strength
						double cMatch = bitMatch(source.cTag, b.cTag);
						
						double distance = dist(source.x, n.x, source.y, n.y);
						if (distanceExponent > 0) {
							distance = std::pow(1 - (distance / 1.4143), distanceExponent); //nomalize to max distance and then square
							cMatch = distance * std::pow(cMatch, cTagExponent);
						}
						else {
							cMatch = std::pow(cMatch, cTagExponent);
						}
						if (cMatch > wireConnectionThreshold) { // make a connection
							//std::cout << "\n  ** connection **   " << nID << "-" << bID << "  to  " << sourceID << std::endl;
							//std::cout << "    distance: " << distance;
							//std::cout << "  cTagMatch: " << bitMatch(source.cTag, b.cTag) << "(" << source.cTag << " " << b.cTag << ")";
							//std::cout << "  cMatch: " << cMatch << std::endl;
							//std::cout << "    wTagMatch: " << bitMatch(source.wTag, b.wTag) << "(" << source.wTag << " " << b.wTag << ")";
							//std::cout << "  final w: " << (bitMatch(source.wTag, b.wTag) * 2.0) - 1.0 << std::endl;
							b.inputConnectionIDs.push_back(sourceID);
							b.inputWeights.push_back((bitMatch(source.wTag, b.wTag)*2.0)-1.0);
						}
					}
				}
			}
		}
	}

	// connect outputs to nearest neurons (possibly, outputs should be branches...)
	void wireOutputs() {
		outputNeuronIDs.clear();
		for (int i = 0; i < nrOutputValues; i++) {
			int outConnection = 0;
			double shortestDistance = dist(outputLocations[i].first, neurons[nrInputValues].x, outputLocations[i].second, neurons[nrInputValues].y);
			for (int sourceID = nrInputValues; sourceID < neurons.size(); sourceID++) {
				const auto& source = neurons[sourceID];
				double distance = dist(outputLocations[i].first, source.x, outputLocations[i].second, source.y);
				if (distance < shortestDistance) {
					shortestDistance = distance;
					outConnection = sourceID;
				}
			}
			outputNeuronIDs.push_back(outConnection);
		}
	}


	///////////////////////////////////////////////////////////////////////////////////
	//output this brain for graphviz //////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	void graphBrain() {

		std::string graphString = "digraph neurons {\n";

		double GV_scale = 10.0;
		double GV_neuronWidth = .05 * GV_scale;
		double GV_neuronHeight = .02 * GV_scale;
		double GV_branchHeight = .02 * GV_scale;
		double GV_fontSize = 10 * GV_scale;

		std::vector<std::string> nNames;
		for (int i = 0; i < neurons.size(); i++) {
			// add input neurons to graph
			

			if (i < nrInputValues) {
				nNames.push_back("in_" + std::to_string(i));
				graphString += nNames.back() + "[width = " + std::to_string(GV_neuronWidth) +
					" height = " + std::to_string(GV_neuronHeight) +
					" pos = \"" + std::to_string(neurons[i].x * GV_scale) + "," + std::to_string(neurons[i].y * GV_scale) +
					"!\" shape = box fontsize = GV_fontSize fillcolor=grey90, style=filled]\n";
			}

			else {
				auto cast_logic = std::dynamic_pointer_cast<Spikey_Neuron_Logic_Device>(neurons[i].logic);
				double dc = cast_logic->deliveryCharge;

				if (cast_logic->defalutFire) {
					// add non-input neurons to graph
					nNames.push_back("n_" + std::to_string(i - (nrInputValues)));

					graphString += nNames.back() + "[label = \"n_"+ std::to_string(i - (nrInputValues)) +
						"__dc_"+std::to_string(dc)+
						"\" width = " + std::to_string(GV_neuronWidth) +
						" height = " + std::to_string(GV_neuronHeight) +
						" pos = \"" + std::to_string(neurons[i].x * GV_scale) + "," + std::to_string(neurons[i].y * GV_scale) +
						"!\" shape = box fontsize = GV_fontSize fillcolor=green, style=filled]\n";
				}
				else {
					nNames.push_back("n_" + std::to_string(i - (nrInputValues)));
					graphString += nNames.back() + "[width = " + std::to_string(GV_neuronWidth) +
						" height = " + std::to_string(GV_neuronHeight) +
						" pos = \"" + std::to_string(neurons[i].x * GV_scale) + "," + std::to_string(neurons[i].y * GV_scale) +
						"!\" shape = box fontsize = GV_fontSize fillcolor=red, style=filled]\n";
				}
			}

		}

		for (int i = nrInputValues; i < neurons.size(); i++) {

			// for each non input neuron...
			std::string nName = nNames[i];
			double numBranch = neurons[i].branches.size();
			double GV_branchWidth = GV_neuronWidth / numBranch;

			for (int j = 0; j < neurons[i].branches.size(); j++) {
				// make a name for each branch
				std::string bName = "b_" + std::to_string(i - (nrInputValues)) + "_" + std::to_string(j);
				if (neurons[i].branches[j].inputConnectionIDs.size() > 0) {
					// if a branch has connections, add to graph
					graphString += bName + "[width = " + std::to_string(GV_branchWidth / 4) + 
						" height = " + std::to_string(GV_branchHeight / 4) +
						" label=\"\" pos=\"" +
						std::to_string((neurons[i].x * GV_scale) + (GV_branchWidth * (j - (numBranch / 2))) + (GV_branchWidth * .5)) + "," +
						std::to_string((neurons[i].y * GV_scale) - (GV_branchHeight * 2)) + "!\" shape = box ]\n";
					// add conneciton from branch to neuron
					graphString += bName + " -> " + nName + " [color=blue]\n";
				}
				for (int k = 0; k < neurons[i].branches[j].inputConnectionIDs.size(); k++) {

					if (neurons[i].branches[j].inputWeights[k] != 0) {
						std::string connectionColor = neurons[i].branches[j].inputWeights[k] >= 0 ? "black" : "red";
						graphString += nNames[neurons[i].branches[j].inputConnectionIDs[k]] + " -> " + bName +
							" [label=" + std::to_string(neurons[i].branches[j].inputWeights[k]) + " fontsize = GV_fontSize penwidth =" +
							std::to_string(std::abs(neurons[i].branches[j].inputWeights[k] * 4) + 1) + " color=" + connectionColor + "]\n";
					}
				}
			}
		}

		for (int i = 0; i < outputNeuronIDs.size(); i++) {
			// make outputs
			std::string oName = "out_" + std::to_string(i);
			std::string nName = nNames[outputNeuronIDs[i]];

			graphString += oName + "[width = " + std::to_string(GV_neuronWidth) + " height = " + std::to_string(GV_neuronHeight) +
				" pos = \"" + std::to_string(outputLocations[i].first * GV_scale) + "," + std::to_string(outputLocations[i].second * GV_scale) +
				"!\" shape = box fontsize = GV_fontSize fillcolor=white, style=filled]\n";

			graphString += nName + "->" + oName + "\n";
		}
		graphString += "}\n";
		FileManager::writeToFile("SWNB.dot",graphString);
	}

	void showBrain() {
		std::stringstream converter; // used to convert strings to other types

		for (int i = 0; i < outputLocations.size(); i++) {
			std::cout << "\noutput" << i << ": " << outputLocations[i].first << " , " << outputLocations[i].second << std::endl;
		}
		for (int i = 0; i < neurons.size(); i++) {
			std::cout << "\nneuron " << i << ": ";
			if (neurons[i].isInputNeuron) {
				std::cout << "input";
			}
			std::cout << "  location: " << neurons[i].x << " , " << neurons[i].y;

			auto cast_logic = std::dynamic_pointer_cast<Spikey_Neuron_Logic_Device>(neurons[i].logic);
			double dc = cast_logic->deliveryCharge;

			std::cout << " with cTag: " << neurons[i].cTag;
			std::cout << "  wTag: " << neurons[i].wTag << std::endl;
			if (!neurons[i].isInputNeuron) {
				std::cout << " logic has  ::  threshold: " << cast_logic->thresholdInit << "  init charge: " << cast_logic->currentChargeInit <<
					"  decay : " << cast_logic->decayRateInit << "  delivery charge: " << cast_logic->deliveryChargeInit << std::endl;
			}
			if (!neurons[i].isInputNeuron) {
				std::cout << "\b   branch count: " << neurons[i].branches.size() << std::endl;
				for (int j = 0; j < neurons[i].branches.size(); j++) {
					std::cout << "     branch: " << j << " has weight: " << neurons[i].branchWeights[j];
					std::cout << " with logic: " << neurons[i].branches[j].logic->serialize().str() << std::endl;
					std::cout << "       cTag: " << neurons[i].branches[j].cTag << "  wTag: " << neurons[i].branches[j].wTag << std::endl;
					for (int k = 0; k < neurons[i].branches[j].inputConnectionIDs.size(); k++) {
						std::cout << "         connection to neuron " << neurons[i].branches[j].inputConnectionIDs[k] << " has weight: " << neurons[i].branches[j].inputWeights[k] << std::endl;
					}
				}
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////////
	// Parts //////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	SWNBrain() = delete;

	SWNBrain(int _nrInNodes, int _nrOutNodes, std::shared_ptr<ParametersTable> PT_ = Parameters::root);

	virtual ~SWNBrain() = default;
	
	virtual void update() override;

	virtual std::shared_ptr<AbstractBrain> makeBrain(std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes) override;
	virtual std::shared_ptr<AbstractBrain> makeCopy(std::shared_ptr<ParametersTable> PT_ = nullptr) override;

	virtual std::shared_ptr<AbstractBrain> makeBrainFrom(std::shared_ptr<AbstractBrain> parent,
		std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes) override;

	virtual std::shared_ptr<AbstractBrain> makeBrainFromMany(std::vector<std::shared_ptr<AbstractBrain>> parents,
		std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes) override;

	virtual std::string description() override;
	virtual DataMap getStats(std::string& prefix) override;
	virtual std::string getType() override { return "Chris"; }

	virtual void resetBrain() override;

	virtual void initializeGenomes(std::unordered_map<std::string, std::shared_ptr<AbstractGenome>>& _genomes);

	virtual std::unordered_set<std::string> requiredGenomes() override {
		return {};
	}

	// convert a brain into data map with data that can be saved to file so this brain can be reconstructed
	// 'name' here contains the prefix that must be so that deserialize can identify relavent data
	virtual DataMap serialize(std::string& name) override;

	// given an unordered_map<string, string> of org data and PT, load data into this brain
	// 'name' here contains the prefix that was used when data was being saved
	virtual void deserialize(std::shared_ptr<ParametersTable> PT, std::unordered_map<std::string, std::string>& orgData, std::string& name) override;

};

inline std::shared_ptr<AbstractBrain> ChrisBrain_brainFactory(
	int ins, int outs, std::shared_ptr<ParametersTable> PT = Parameters::root) {
	return std::make_shared<SWNBrain>(ins, outs, PT);
}
