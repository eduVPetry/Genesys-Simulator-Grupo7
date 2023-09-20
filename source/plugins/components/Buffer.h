/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Buffer.h
 * Author: rafael.luiz.cancian
 *
 * Created on
 */

#ifndef BUFFER_H
#define BUFFER_H

#include "../../kernel/simulator/ModelComponent.h"
#include "../data/Queue.h"
#include "../data/SignalData.h"

/*!
 This component ...
 */
class Buffer : public ModelComponent {
public:
	enum class AdvanceOn : int {
		NewArrivals = 0, Time = 1, Signal = 2
	};
	enum class ArrivalOnFullBufferRule : int {
		Dispose = 0, SendToBulkPort = 1, ReplaceLastPosition = 2
	};
public: // constructors
	Buffer(Model* model, std::string name = "");
	virtual ~Buffer() = default;
public: // virtual
	virtual std::string show();
public: // static
	static PluginInformation* GetPluginInformation();
	static ModelComponent* LoadInstance(Model* model, PersistenceRecord *fields);
	static ModelDataDefinition* NewInstance(Model* model, std::string name = "");
public:
	Buffer::ArrivalOnFullBufferRule getarrivalOnFullBufferRule() const;
	void setArrivalOnFullBufferRule(Buffer::ArrivalOnFullBufferRule newArrivalOnFullBufferRule);
	Buffer::AdvanceOn getadvanceOn() const;
	void setAdvanceOn(Buffer::AdvanceOn newAdvanceOn);
	unsigned int getcapacity() const;
	void setCapacity(unsigned int newCapacity);
	Queue *getinternalQueue() const;
	void setInternalQueue(Queue *newInternalQueue);
	SignalData *getsignal() const;
	void setSignal(SignalData *newSignal);

	std::string getadvanceTimeExpression() const;
	void setAdvanceTimeExpression(const std::string &newAdvanceTimeExpression, Util::TimeUnit timeunit = Util::TimeUnit::unknown);
	Util::TimeUnit getadvanceTimeTimeUnit() const;
	void setAdvanceTimeTimeUnit(Util::TimeUnit newAdvanceTimeTimeUnit);

protected: // must be overriden
	virtual bool _loadInstance(PersistenceRecord *fields);
	virtual void _saveInstance(PersistenceRecord *fields, bool saveDefaultValues);
	virtual void _onDispatchEvent(Entity* entity, unsigned int inputPortNumber);
protected: // could be overriden by derived classes
	virtual bool _check(std::string* errorMessage);
	/*! This method returns all changes in the parser that are needed by plugins of this ModelDatas. When connecting a new plugin, ParserChangesInformation are used to change parser source code, whch is after compiled and dinamically linked to to simulator kernel to reflect the changes */
	virtual ParserChangesInformation* _getParserChangesInformation();
	virtual void _initBetweenReplications();
	/*! This method is necessary only for those components that instantiate internal elements that must exist before simulation starts and even before model checking. That's the case of components that have internal StatisticsCollectors, since others components may refer to them as expressions (as in "TVAG(ThisCSTAT)") and therefore the modeldatum must exist before checking such expression */
	virtual void _createInternalAndAttachedData(); /*< A ModelDataDefinition or ModelComponent that includes (internal) ou refers to (attach) other ModelDataDefinition must register them inside this method. */
	virtual void _addProperty(PropertyBase* property);
private: // methods
	unsigned int _handlerForSignalDataEvent(SignalData* signalData);
private: // attributes 1:1

	const struct DEFAULT_VALUES {
		const AdvanceOn advanceOn = AdvanceOn::NewArrivals;
		const ArrivalOnFullBufferRule arrivalOnFullBufferRule = ArrivalOnFullBufferRule::Dispose;
		const unsigned int capacity = 1;
		const std::string advanceTimeExpression = "1";
		const Util::TimeUnit advanceTimeTimeUnit = Util::TimeUnit::second;
	} DEFAULT;
	ArrivalOnFullBufferRule _arrivalOnFullBufferRule = DEFAULT.arrivalOnFullBufferRule;
	AdvanceOn _advanceOn = DEFAULT.advanceOn;
	unsigned int _capacity = DEFAULT.capacity;
	std::string _advanceTimeExpression = DEFAULT.advanceTimeExpression;
	Util::TimeUnit _advanceTimeTimeUnit = DEFAULT.advanceTimeTimeUnit;
	Queue* _internalQueue = nullptr;
	SignalData* _attachedSignal = nullptr;
private: // attributes 1:n
};

#endif /* BUFFER_H */
