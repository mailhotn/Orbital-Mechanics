# Orbital-Mechanics

This repository includes a number of classes and functions used to simulate and evaluate the performance of constellations of satellies.
The main focus of this code base is constellations for TDOA based geolocation.

## The Constellation class & subclasses
In order to simulate a constellation's dynamics, the constellation must be defined using one of the available constellation classes.
Each such class is a subclass of the Constellation superclass. This superclass contains the properties that are universal to all types of constellation.
Each Constellation class must have a certain set of methods which create the initial state of the constellation.
Each Constellation class must have the following methods:
- _InitialOeOsc_
- _InitialStateEci_
- _InitialOeMean_

## The Propagator class
In addition to the Constellation, a Propagator must be created using the appropriate class.
The dynamics of the constellation can then be simulated over a given time period using the "Prop" methods:
- _PropEciTb_ propagates the state of the satellites according to Keplerian Mechanics.
- _PropEciJ2_ propagates the state of the satellites with the J<sub>2</sub> perturbation.
- _PropOeMean_ propagates the mean orbital elements of the constellation using a numerical solver, with the J<sub>2</sub> perturbation.
- _PropOeMeanFast_ propagates the mean orbital elements of the constellation without using a numerical solver, with the J<sub>2</sub> perturbation.
- _PropOeOsc_ propagates the osculating orbital elements of the constellation with the J<sub>2</sub> perturbation.