

function simulation( C, Cim, Cs, conf, Ct, Gi, stats ){

	this.C = C
	this.Cim = Cim
	this.Cs = Cs
	this.conf = conf
	this.maxTCells = this.conf["MAX_TCELLS"]
	this.Gi = Gi
	this.stats = stats
	this.math = require("mathjs")

	// max x,y,z coordinates in the grid
	this.maxc = [this.C.field_size.x-1, this.C.field_size.y-1 ]
	if( this.C.ndim == 3 ){
		this.maxc[2] = this.C.field_size.z-1
	}

	this.viewtracks = this.conf["VIEWTRACKS"]
	if( this.viewtracks ){
		this.Ct = Ct
	}

	this.save = this.conf["SAVEIMG"]
	this.runtime=this.conf["RUNTIME"]
	this.maxKilling=this.conf["MAXKILLING"]
	this.simtype = this.conf["SIMTYPE"] || "2D"
	this.entryBias = this.conf["ENTRYBIAS"]
	this.entryBiasStrength = 1

	// to track the time of the actual simulation (no burnin)
	this.time = 0
	this.borderthreshold = 30
	this.Tcell_infection_borders = []
	this.infectionChance = 0.000025
	this.avg_border_CTL_infection = 0
	this.borderingparameter = 2

	this.infectionlist = []
}

simulation.prototype = {

	// Seed cells and allow a burnin phase, then start the simulation.
	initialize : function(){

		var nrcells = this.conf["NRCELLS"], burnin = this.conf["BURNIN"]
		var cellkind, i, stromavoxels

		// Seed the right number of cells for each cellkind
		var totalcells = 0
		for( cellkind = 0; cellkind < nrcells.length; cellkind ++ ){

			for( i = 0; i < nrcells[cellkind]; i++ ){
				let newid = this.Gi.seedCell( cellkind+1 )
				if(cellkind == 3){
					this.infectionlist[newid] = 1
				}
				totalcells++
			}
		}

		// Simulate the burnin phase
		for( i = 0; i < burnin; i++ ){
			this.C.monteCarloStep()
		}

		// Call timestep to start the real simulation.
		// this.timestep()

	},

	// returns cells of type "type" that border "affectors"
	getCellsToAffect : function( affectors, type ){
		var cellsToAffect = {}
		let cbpig = this.C.cellBorderPixelIndices()
		let cbpi = {}
		for(let i = 0; i < this.C.t2k.length; i++){
			cbpi[i] = []
		}
		for(let t of cbpig){
			let n = cbpig.next().value
			if(n == undefined){
				break
			}
			cbpi[n[1]].push(n[0])
		}
		for ( i = 0; i < affectors.length; i++ ) {
			infectedCellNeighbors = this.stats.cellNeighborsList(affectors[i], cbpi)[1]
			// console.log(this.stats.cellNeighborsList(affectors[i], cbpi))
			for (const [key, value] of Object.entries(infectedCellNeighbors)) {
				if (this.C.t2k[key] == type){
					if ( key in cellsToAffect ){
						cellsToAffect[key] = cellsToAffect[key] + value
					}
					else {
						cellsToAffect[key] = value
					}
				}
			}
		}
		return cellsToAffect
	},

	// returns cells of type "type"
	findCells : function(type){
		var cells = []
		for( i = 0; i < this.C.t2k.length; i++ ) {
			if( this.C.t2k[i] == type ) {
				cells.push(i)
			}
		}
		return cells
	},

	// CTLs neighboring infected cells have a chance to kill them
	kill : function(){
		var killers = this.findCells(1).concat(this.findCells(5))
		var infectedNeighbors = this.getCellsToAffect( killers, 3 )

		for (const [key, value] of Object.entries(infectedNeighbors)){
			let rand = Math.random()
			if (rand < value * (1/this.maxKilling) ){
				this.C.setCellKind(key, 4)
			}
		}
	},

	// infected cells neighboring healthy skin cells have a chance ot infect them
	infectOthers : function(){
		var infected = this.findCells(3)
		var skinNeighbors = this.getCellsToAffect( infected, 2 )

		for (const [key, value] of Object.entries(skinNeighbors)){
			let rand = Math.random()
			if (rand < value * this.infectionChance){
				this.infectionlist[key] = 1
				this.C.setCellKind(key, 3)
			}
		}
	},

	// infected skin cells will get more infected over time (this only has a consequence for their color)
	getMoreInfected : function(){
		for( i = 0; i < this.C.t2k.length; i++ ) {
			if( this.infectionlist[i] > 0 && this.infectionlist[i] < 2000 ) {
				this.infectionlist[i] = this.infectionlist[i] + 1
			}
		}
	},

	// if a CTL replaces an infected cell, a healthy skin cell bordering the infection should be infected immediately
	replaceInfectionBorderCell : function(cell_idToReplace){
		// replace border cell with most infected neighbors with this cell
		var skinNeighbors = this.getCellsToAffect( this.findCells(3), 2 )
		var maxKey = -1
		var maxVal = -1
		for (const [key, value] of Object.entries(skinNeighbors)){
			if (value > maxVal) {
				maxKey = key
				maxVal = value
			}
		}
		this.infectionlist[maxKey] = this.infectionlist[cell_idToReplace]
		this.C.setCellKind(maxKey, 3)
	},

	nmod : function(x, N) {
		return ((x % N) + N) % N;
	},

	t21 : function(x,y,N){
		return this.nmod(y,N)*N+this.nmod(x,N)
	},

	// return x and y based on chemokine level (biased entry)
	biasedEntry : function() {
		let xpos = Math.floor( Math.random()*( this.C.field_size.x ) )
		let ypos = Math.floor( Math.random()*( this.C.field_size.y ) )
		let maxConcentration = this.math.max(this.C.chemokinelevel)
		let sumConcentrations =  this.math.sum(this.C.chemokinelevel)
	  let maxProb = maxConcentration/sumConcentrations
		
		let thisProb = this.C.chemokinelevel.get([this.t21(this.math.floor(xpos),this.math.floor(ypos),(this.C.field_size.x)),0])/sumConcentrations
		var thisRealProb = (thisProb/maxProb)*this.entryBiasStrength + (1-this.entryBiasStrength)

		while ( Math.random() > thisRealProb ) {
			xpos = Math.floor( Math.random()*( this.C.field_size.x ) )
			ypos = Math.floor( Math.random()*( this.C.field_size.y ) )
			thisProb = this.C.chemokinelevel.get([this.t21(this.math.floor(xpos),this.math.floor(ypos),(this.C.field_size.x)),0])/sumConcentrations
			thisRealProb = (thisProb/maxProb)*this.entryBiasStrength + (1-this.entryBiasStrength)
		}
		return [xpos, ypos]
	},

	// adds CTL to the grid
	addTCell : function() {
		var cell_idToReplace = -1
		typeToReplace = 1
		while ( typeToReplace != 2 && typeToReplace != 3 ) {
			let x = Math.floor( Math.random()*( this.C.field_size.x ) )
			let y = Math.floor( Math.random()*( this.C.field_size.y ) )
			if ( this.entryBias == 1 ) {
				let xy = this.biasedEntry()
				x = xy[0]
				y = xy[1]
			}
			cell_idToReplace = this.C.pixt([x, y])
			typeToReplace = this.C.cellKind(cell_idToReplace)
		}
		if (this.infectionlist[cell_idToReplace] > 0) {
			this.replaceInfectionBorderCell(cell_idToReplace)
		}
		this.infectionlist[cell_idToReplace] = -1
		this.C.setCellKind(cell_idToReplace, 1)
	},

	// changes or resets CTL behavior based on whether or not they border infected cells enough
	changeCTLBehavior : function() {
		let CTLs = this.findCells(1).concat(this.findCells(5))
		// let cbpi = this.Cs.cellborderpixelsi()
		let cbpig = this.C.cellBorderPixelIndices()
		let cbpi = {}
		for(let i = 0; i < this.C.t2k.length; i++){
			cbpi[i] = []
		}
		for(let t of cbpig){
			let n = cbpig.next().value
			if(n == undefined){
				break
			}
			cbpi[n[1]].push(n[0])
		}
		for(let i = 0; i < CTLs.length; i++) {
			let neighbors = this.stats.cellNeighborsList(CTLs[i], cbpi)
			let infectedNeighbor = false
			let border = 0
			for(let j = 0; j < neighbors.length; j++) {
				if(this.C.cellKind(neighbors[0][j]) == 3 ) {
					border += neighbors[1][neighbors[0][j]]
				}
			}
			if(border > this.avg_border_CTL_infection*(this.borderingparameter/3)) {
				infectedNeighbor = true
			}
			if(this.C.cellKind(CTLs[i]) == 1 && infectedNeighbor) {
				this.C.setCellKind(CTLs[i], 5)
			}
			if(this.C.cellKind(CTLs[i]) == 5 && !infectedNeighbor) {
				this.C.setCellKind(CTLs[i], 1)
			}
		}
	},

	// Update the grid and either draw it or print track output
	timestep : function(){
		// CTLs first might kill infected cells
		this.kill()
		// infected cells now may affect others
		this.infectOthers()
		// infected cells will get more infected
		this.getMoreInfected()
		// sometimes a CTL will appear
		if ((this.time + 1) % 48 == 0 && this.Cs.countCells(1) + this.Cs.countCells(5) < this.maxTCells) {
			this.addTCell()
		}
		// performs a MCS
		this.C.monteCarloStep()
		// chenges behavior of CTLs if needed
		this.changeCTLBehavior()
		this.time++
	},

	atBorder : function(){

		var cbp = this.C.cellborderpixels
		var i,p,j

		// only check in detail if centroid is within 20px of grid borders
		if( this.nearBorder(this.borderthreshold) ){

			// loop over the cellborderpixels, and check if they are at the grid border.
			for( i = 0; i < cbp.length; i++ ){
				p = this.C.i2p( cbp.elements[i] )
				for( j = 0; j < p.length; j++ ){
					if( p[j] == 0 || p[j] == this.maxc[j] ){
						return true
					}
				}
			}

		}
		return false

	},

	nearBorder : function( threshold ){

		// since there is only one cell, this is an array
		// of length 1 containing an object with "id", x,y,z
		// of the cell t = "id".
		// Still, loop for consistency:
		var centroid = this.Cs.getCentroids()
		var i,c


		for( i = 0; i < centroid.length; i++ ){
			c = centroid[i]
			var xmax = this.C.field_size.x - threshold - 1,
				ymax = this.C.field_size.y - threshold -1,
				zmax = this.C.field_size.z-threshold-1

			if( c.x < threshold || c.x > xmax ){
				return true;
			}
			if( this.simtype != "1D" && ( c.y < threshold || c.y > ymax ) ){
				return true;
			}
			if( this.C.ndim == 3 ){
				if( c.z < threshold || c.z > zmax ){
					return true;
				}
			}
		}

		return false
	},


	// Draw the grid
	drawCanvas : function(){

		var canvascolor=this.conf["CANVASCOLOR"], stromacolor=this.conf["STROMACOLOR"],
			showborders=this.conf["SHOWBORDERS"],
			cellcolor=this.conf["CELLCOLOR"], actcolor=this.conf["ACTCOLOR"],
			trackcolor=this.conf["TRACKCOLOR"],
			nrcells=this.conf["NRCELLS"],cellkind

		// Clear canvas and draw stroma border
		this.Cim.clear( canvascolor )

		// Draw each cellkind appropriately
		for( cellkind = 0; cellkind < nrcells.length; cellkind ++ ){
			if( cellcolor[ cellkind ] != -1 ){
				this.Cim.drawCells( cellkind+1, cellcolor[cellkind])//, this.infectionlist )
			}
			if( actcolor[ cellkind ] ){
				this.Cim.drawActivityValues( cellkind + 1 )
			}
			if( showborders[ cellkind ] ){
				//this.Cim.drawCellBorders( cellkind+1, "AAAAAA" )
				this.Cim.drawCellBorders( cellkind+1, "000000" )
			}
		}
	},

	reset : function(){


		// remove all cells from the grid, but keep stroma
		var cellpix1 = this.C.cellpixelstype, cellpix = Object.keys( cellpix1 ), i
		for( i = 0; i < cellpix.length; i++ ){
			if( cellpix1[i] != 0 ){
				this.C.delpixi( cellpix[i] )
			}
		}
		// Also remove the cellborderpixels
		this.C.cellborderpixels = new DiceSet()
		// call initialize to re-seed cells
		this.initialize()
	}


}


/* This allows using the code in either the browser or with nodejs. */
if( typeof module !== "undefined" ){
	// var DiceSet = require("../../cpm/src/DiceSet.js")
	module.exports = simulation
}
