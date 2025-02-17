/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   StatisticsDefaultImpl1.cpp
 * Author: rafael.luiz.cancian
 *
 * Created on 1 de Agosto de 2018, 21:03
 */

#include <complex>

#include "StatisticsDefaultImpl1.h"
#include "../TraitsKernel.h"
//#include "Integrator_if.h"
//#include "ProbDistribDefaultImpl1.h"

//using namespace GenesysKernel;

StatisticsDefaultImpl1::StatisticsDefaultImpl1() {
	//_collector = new TraitsKernel<Statistics_if>::CollectorImplementation();
	_collector = new TraitsKernel<Model>::StatisticsCollector_CollectorImplementation();
	_collector->setAddValueHandler(setCollectorAddValueHandler(&StatisticsDefaultImpl1::collectorAddHandler, this));
	_collector->setClearHandler(setCollectorClearHandler(&StatisticsDefaultImpl1::collectorClearHandler, this));
	//_collector->setAddValueHandler(std::bind(&StatisticsDefaultImpl1::collectorAddHandler, this, std::placeholders::_1));
	this->initStatistics();
}

StatisticsDefaultImpl1::StatisticsDefaultImpl1(Collector_if* collector) {
	_collector = collector;
	_collector->setAddValueHandler(setCollectorAddValueHandler(&StatisticsDefaultImpl1::collectorAddHandler, this));
	_collector->setClearHandler(setCollectorClearHandler(&StatisticsDefaultImpl1::collectorClearHandler, this));
	//_collector->setAddValueHandler(std::bind(&StatisticsDefaultImpl1::collectorAddHandler, this, std::placeholders::_1));
	this->initStatistics();
}

void StatisticsDefaultImpl1::collectorAddHandler(double newValue, double newWeight) {
	_elems = _collector->numElements();
	if (newValue < _min) {
		_min = newValue;
	}
	if (newValue > _max) {
		_max = newValue;
	}

	// alternative 1
	// equally
	//_sumData += newValue;
	//double oldAverage = _average;
	//_average = _sumData / _elems;
	//_sumDataSquare += (newValue-oldAverage)*(newValue-_average);
	//_unweightedvariance = _sumDataSquare/_elems;
	//if (_elems>1) {
	//	_unbiasedVariance =_sumDataSquare/(_elems-1);
	//}
	// Considering weight
	_sumWeight += newWeight;
	_sumWeightSquare += newWeight*newWeight;
	double oldAverage = _average;
	_average = oldAverage + (newWeight/_sumWeight)*(newValue-oldAverage);
	_sumData += newWeight*(newValue-oldAverage)*(newValue-_average);
	_variance = _sumWeight==1 ? 0 :_sumData/(_sumWeight-1);
	_unbiasedVariance = _sumWeight==0 ? 0 : _sumData / (_sumWeight-_sumWeightSquare/_sumWeight);

	// alternative 2 (numerical instability)
	//_average = _average + (newValue - _average)/_elems;  // or bellow
	//_average = (_average * (elems - 1) + newValue) / elems;  // this approach propagates the numeric error
	//_variance = (_variance * (elems - 1) + pow(newValue - _average, 2)) / elems;  // this approach propagates the numeric error

	_stddeviation = sqrt(_variance);
	_variationCoef = (_average != 0 ? _stddeviation / _average : 0.0);
	_halfWidth = _criticalTn_1 * (_stddeviation / std::sqrt(_elems));
}

void StatisticsDefaultImpl1::collectorClearHandler() {
	this->initStatistics();
}

void StatisticsDefaultImpl1::initStatistics() {
	_elems = 0;
	_min = +1e+99; //@TODO: Change by the double min and man
	_max = -1e+99;
	_sumData = 0.0;
	_sumWeight = 0.0;
	_sumWeightSquare = 0.0;
	_sumDataSquare = 0.0;
	_average = 0.0;
	_variance = 0.0;
	_unweightedvariance = 0.0;
	_unbiasedVariance = 0.0;
	_stddeviation = 0.0;
	_variationCoef = 0.0;
	_halfWidth = 0.0;
}

unsigned int StatisticsDefaultImpl1::numElements() {
	return this->getCollector()->numElements();
}

double StatisticsDefaultImpl1::min() {
	if (_elems > 0)
		return _min;
	else
		return 0.0;
}

double StatisticsDefaultImpl1::max() {
	if (_elems > 0)
		return _max;
	else
		return 0.0;
}

double StatisticsDefaultImpl1::average() {
	return _average;
}

double StatisticsDefaultImpl1::variance() {
	return _variance;
}

double StatisticsDefaultImpl1::stddeviation() {
	return _stddeviation;
}

double StatisticsDefaultImpl1::variationCoef() {
	return _variationCoef;
}

double StatisticsDefaultImpl1::halfWidthConfidenceInterval() {
	return _halfWidth;
}

void StatisticsDefaultImpl1::setConfidenceLevel(double confidencelevel) {
	_confidenceLevel = confidencelevel;
	//Integrator_if* integrator = new TraitsKernel<Integrator_if>::Implementation();
	_criticalTn_1 = 1.96; //integrator->integrate()

}

double StatisticsDefaultImpl1::confidenceLevel() {
	return _confidenceLevel;
}

unsigned int StatisticsDefaultImpl1::newSampleSize(double halfWidth) {
	return 0;
}

Collector_if* StatisticsDefaultImpl1::getCollector() const {
	return this->_collector;
}

void StatisticsDefaultImpl1::setCollector(Collector_if* collector) {
	this->_collector = collector;
}
