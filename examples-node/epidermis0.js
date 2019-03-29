// require("../src/DiceSet.js")
// require("../src/CPM.js")
// var CPMStats = require("../src/CPMStats.js")
// var CPMCanvas = require("../src/CPMCanvas.js")
// var BGCanvas = require("../src/BGCanvas.js")
// var CPMchemotaxis = require("../src/CPMchemotaxis.js")
// var TrackCanvas = require("../src/TrackCanvas.js")
var simulation = require("./simulation-tissue.js")
var simsettings = {
	SIMTYPE : "2D",
	NRCELLS : [0,280,0,0,0],
	BURNIN : 50,
	RUNTIME : 12000,
	CANVASCOLOR : "FFFFFF",
	STROMACOLOR : "AAAAAA",
	CELLCOLOR : ["00CC00","4A4C4F","9A2222","222222","77CC77"],
	ACTCOLOR : [false,false,false,false,false],
	VIEWTRACKS : false,
	TRACKCOLOR : ["FF0000"],
	SHOWBORDERS : [false,true,true,true,true],
	SAVEIMG : false,
	FRAMERATE : 10,
	SAVEPATH : "movies",
	MCSRATE : 5,
	MAX_TCELLS : 100,
	MAXKILLING : 1,
	ENTRYBIAS : 0
}

// var CsetChemotaxis = require("./CPMskin-template.js")
// var simsettings = require("./CPMskin-settings.js")
var CPM = require("../build/cpm-cjs.js")
const fs = require('fs')
var runtime = parseInt(process.argv[2]) || 1000 /////////////////////////////////////////
var savetime = parseInt(process.argv[3]) || 10 //////////////////////////////////////////
var field_size = parseInt(process.argv[4]) || 200 //////////////////////////////////////
var kera_cells_factor = (field_size*field_size)/(200*200) //////////////////////////////
simsettings["MAX_TCELLS"] = 100*( (field_size*field_size) / (600*600) ) ////////////////
var chemotaxis = parseInt(process.argv[5]) || 0 ////////////////////////////////////////
var infectionStart = parseFloat(process.argv[6]) || 0.1 ////////////////////////////////
var entryBias = parseInt(process.argv[7]) || 0 /////////////////////////////////////////
simsettings["ENTRYBIAS"] = entryBias ///////////////////////////////////////////////////
var killingTime = parseInt(process.argv[8]) || 15 //15, 30, 45, 60 /////////////////////
var avg_border_CTL_infection = 30.46 ///////////////////////////////////////////////////
var simulationType = parseInt(process.argv[9]) || 1 ////////////////////////////////////
var infectionChance = parseFloat(process.argv[10]) || 0.000025 /////////////////////////
var maxAct = parseInt(process.argv[11]) || 20 //////////////////////////////////////////
var lAct = parseInt(process.argv[12]) || 1000 //////////////////////////////////////////
var borderingparameter = parseInt(process.argv[13]) || 2 ///////////////////////////////
var idnr = parseInt(process.argv[14]) || 0 /////////////////////////////////////////////
var report = true //////////////////////////////////////////////////////////////////////

var C,Gi,Cim,Ct,Cs, stopped=true, zoom=4, wrap=[0,0], stats

function initialize(){
	let C = new CPM.CPM( [field_size,field_size], {
		T : 20,
		torus : true
	})
	Cs = new CPM.PostMCSStats()
	stats = new CPM.Stats(C)
	C.add( Cs )
	C.add( new CPM.Adhesion( { J: [[0,0,0,0,0,0], [0,100,20,20,20,100], [20,20,40,40,20,20],
		[20,20,20,40,20,20], [20,20,50,50,50,20], [0,100,20,20,20,100]] } ) )
	C.add( new CPM.VolumeConstraint( { V : [0,100,145,145,145,100], LAMBDA_V : [0,50,50,50,50,50] } ) )
	C.add( new CPM.PerimeterConstraint( { P : [0,125,145,145,145,125], LAMBDA_P : [0,2,2,2,2,2] } ) )

	// Set different kinds of behaviors
	let lambda_chem_1 = chemotaxis
	let lambda_chem_5 = chemotaxis
	let lambda_act_1 = lAct
	let lambda_act_5 = lAct
	let max_act_1 = maxAct
	let max_act_5 = maxAct
	if(simulationType == 2){
		lambda_chem_5 = 0
		lambda_act_5 = 0
	}
	if(simulationType == 3){
		lambda_chem_5 = 0
	}

	C.add( new CPM.ActivityConstraint ( { LAMBDA_ACT : [0,lambda_act_1,0,0,0,lambda_act_5], MAX_ACT : [0,max_act_1,0,0,0,max_act_5], ACT_MEAN : "geometric" } ) )
	let chemokinegradient = new CPM.ChemotaxisConstraint (
		{SECRETOR: 3, RESOLUTION_DECREASE: 10, MM_PER_PIXEL: .38/600, SECOND_PER_MCS: 1, D: 6.2 * Math.pow(10, -5),
			DIFFUSION_PER_MCS: 10, SECRETION: 100, DECAY: .15, LAMBDA_CHEMOTAXIS : [0,lambda_chem_1,0,0,0,lambda_chem_5]
	} )
	C.add( chemokinegradient )
	C.chemokinelevel = chemokinegradient.chemokinereal

	Gi = new CPM.GridInitializer(C)
	Cim = new CPM.Canvas( C, {zoom:zoom} )

	simsettings["MAXKILLING"] = killingTime * 60 * avg_border_CTL_infection
	simsettings["RUNTIME"] = runtime
	simsettings["NRCELLS"][1] *= kera_cells_factor
	sim = new simulation( C, Cim, Cs, simsettings, Ct, Gi, stats )
	sim.initialize()
	sim.borderingparameter = borderingparameter
	sim.infectionChance = infectionChance
	sim.avg_border_CTL_infection = avg_border_CTL_infection
	seedInfection(C)
	sim.drawCanvas()
	logData()
	report = true
	startSim()
}

function seedInfection(C) {
	let centroids = []
	for(let t of C.cellIDs()) {
		centroids.push(Cs.centroidWithTorusCorrection(t))
	}
	for (let i = 0; i < centroids.length; i++) {
		if(Math.sqrt(Math.pow((centroids[i][0] - (field_size/2)), 2) + Math.pow((centroids[i][1] - (field_size/2)), 2)) < (field_size*infectionStart)) {
			sim.infectionlist[i] = 2000 / 2
			C.setCellKind(i, 3)
		}
	}
}

function logData() {
	console.log(sim.time, chemotaxis, killingTime, entryBias, Cs.countCells(1), Cs.countCells(5), Cs.countCells(2), Cs.countCells(3), Cs.countCells(4) )
}

// Perform MCSs in the model
function step(){
	while ( sim.time <= sim.runtime && !sim.stop ) {
		// MCS is performed
		sim.timestep()
		if (savetime != 0) {
			if (sim.time % savetime == 0) {
				sim.drawCanvas()
				Cim.writePNG("output/" + sim.time + "_" + chemotaxis + "_" + killingTime + "_" + simulationType + "_" + idnr + ".png")
			}
		}
		//log data every 30 MCS and stop simulation if no more infected cells are left
		if(report && sim.time % 30 == 0){
			logData()
			if(Cs.countCells(3) == 0){
				stopSim()
			}
		}
	}
}

// For controlling the simulation
function startSim(){
	sim.stop = false
	step()
}

function stopSim(){
	sim.stop = true
}

initialize()
