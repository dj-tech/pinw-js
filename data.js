
//margini
var margin_isoform = {top: 100, right: 15, bottom: 15, left: 15};
//dimensione della finestra di visualizzazione dell'isoforma
var height = window.innerHeight + 100 - margin_isoform.top - margin_isoform.bottom,
width = window.innerWidth - margin_isoform.left - margin_isoform.right;

/* EXONS_STRUCTURE
 * extract_exons -> dati degli esoni estratti dal file json
 * extract_regions -> dati delle regioni estratti dal file json
 * 
 * Crea la struttura dei dati per gli esoni. Concatena le sequenze
 * nucleotidiche delle regioni e aggiunge la scala di colori per 
 * rappresentare i vari esoni. Riporta anche gli ID delle regioni
 * che compongono l'esone.
 */
function exons_structure(extract_exons, extract_regions, extract_boundary){
	
	l = extract_exons.length;
	
	var exons = [];
	
	for(i = 0; i < l; i++){
		
		var reg = [];
		
		//left & right boundary
		l_b = extract_exons[i].left_boundary;
		r_b = extract_exons[i].right_boundary;
		
		var seq = "";
		var flag_seq = false;
		var flag_alt = false;
		
		//regione boundary di sinistra
		region_left = extract_regions[extract_boundary[l_b].first + 1];
		//regione boundary di destra
		region_right = extract_regions[extract_boundary[r_b].first];
			
		//tipologia regione di sinistra	
		type_region_left = region_left.type;
		//tipologia regione di destra
		type_region_right = region_right.type;
		
		//esclude le regioni non codificanti
		if((type_region_left == "codifying") & (type_region_left != "unknow")){
			start_exon = region_left.start;
			if((type_region_right == "codifying") & (type_region_right != "unknow"))
				end_exon = region_right.end;
			
			//assembla la sequenza nucleotidica e salva le regioni appartenenti all'esone
			for(j = region_left.id, k = 0; j <= region_right.id; j++, k++){
			
				if(extract_regions[j].sequence == null)
					flag_seq = true;
				else{
					if(extract_regions[j].type == "codifying")
						seq = seq.concat(extract_regions[j].sequence);
					else 
						flag_seq = true;
					}
			
				reg[k] = extract_regions[j].id;
				if(extract_regions[j].alternative == true)
					flag_alt = true;
			}
			
			if(flag_seq == true)
				seq = null;
				
			
			
			//costruisce l'oggetto esone
			exons.push({
					"id" : i,
					"start" : start_exon,
					"end" : end_exon,
					"sequence" : seq,
					"regions" : reg,
					"alternative" : flag_alt
			});
			reg = [];
		
		}
	}
	
	return exons;	
}

/* INTRONS_STRUCTURE
 * extract_introns -> dati degli introni estratti dal file json
 * extract_regions -> dati delle regioni estratti dal file json
 * 
 * Crea la struttura dei dati per gli introni. Concatena le sequenze
 * nucleotidiche delle regioni e aggiunge la scala di colori per 
 * rappresentare i vari introni. Riporta anche gli ID delle regioni
 * che compongono l'introne. Calcola il pattern dai suffissi e prefissi
 * delle regioni che compongono l'introne.
 */
function introns_structure(extract_introns, extract_regions, extract_boundary){
	
	l = extract_introns.length;
	var introns = [];
	var reg = [];

	
	for(i = 0; i < l; i++){
		
		var reg = [];
		var seq = "";
		var flag_seq = false;
		flag_intron_ok = true;
		
		//left & right boundary
		l_b = extract_introns[i].left_boundary;
		r_b = extract_introns[i].right_boundary;
		
		
		
		//regione boundary di sinistra
		region_left = extract_regions[extract_boundary[l_b].first + 1];
		//regione boundary di destra
		region_right = extract_regions[extract_boundary[r_b].first];
			
		//tipologia regione di sinistra	
		type_region_left = region_left.type;
		alternative_region_left = region_left.alternative;
		//tipologia regione di destra
		type_region_right = region_right.type;
		alternative_region_right = region_right.alternative;
		
		//esclude le regioni non codificanti
		if((type_region_left != "unknow") & (type_region_right != "unknow")){
			if((extract_boundary[l_b].type ==  "5") | (extract_boundary[l_b].type == "both")){
				if((type_region_left == "codifying") & (alternative_region_left == true))
					start_intron = region_left.start;
				else
					if(type_region_left == "spliced")
						start_intron = region_left.start;
					else
						flag_intron_ok == false;
			}
			else
				flag_intron_ok = false;
			
			if((extract_boundary[r_b].type ==  "3") | (extract_boundary[r_b].type == "both")){
				if((type_region_right == "codifying") & (alternative_region_right == true))
					end_intron = region_right.start;
				else
					if(type_region_right == "spliced")
						end_intron = region_right.end;
					else
						flag_intron_ok = false;
			}
			else
				flag_intron_ok = false;
		}
		else
			flag_intron_ok = false;
			
		if(flag_intron_ok){
			//assembla la sequenza nucleotidica e salva le regioni appartenenti all'introne
			for(j = region_left.id, k = 0; j <= region_right.id; j++, k++){
			
				if(extract_regions[j].sequence == null)
					flag_seq = true;
				else{
					if(((extract_regions[j].type == "codifying") & (extract_regions[j].alternative == true)) | extract_regions[j].type == "spliced")
						seq = seq.concat(extract_regions[j].sequence);
					else
						flag_seq = true;
					}
			
				reg[k] = extract_regions[j].id;
			}
			
			if(flag_seq == true)
				seq = null;
		
			l_suffix = extract_introns[i].suffix.length;
			pattern = extract_introns[i].prefix.substr(0, 2).concat(extract_introns[i].suffix.substr(l_suffix - 2, l_suffix));
			
			introns.push({
					"start" : start_intron,
					"end" : end_intron,
					"sequence" : seq,
					"pattern" : pattern,
					"regions" : reg
			});
			reg = [];
			
		}
		
	}
	
	return introns;	
}

/* SPLICE_SITES_STRUCTURE
 * extract_boundaries -> dati dei boundary estratti dal file json
 * extract_regions -> dati delle regioni estratti dal file json
 * 
 * Crea la struttura dei dati per gli splice_sites. Riporta la posizione
 * e la tipologia di ogni splice_sites 
 */
function splice_site_structure(extract_boundaries, extract_regions){
	
	l = extract_boundaries.length;
	var s_s = [];
	
	for(i = 0; i < l; i++){
		
		if(extract_boundaries[i].first == -1){
			p = null;
			t = "unknow";
		}
			
		else{
			t = extract_boundaries[i].type;
			
			if((t == "5") | (t == "both"))
				p = extract_regions[extract_boundaries[i].first + 1].start;
			if((t == "3") | (t == "term"))
				p = extract_regions[extract_boundaries[i].first].end;
			if(t == "init")
				p = extract_regions[extract_boundaries[i].first].start;
				
			s_s.push({
				"position" : p,
				"type" : t			
			});
		}
	}
	
	return s_s;
}

/* ISOFORM_SCALE
 * e_r -> dati delle regioni estratti dal file json
 * 
 * Definisce un range per trasformare la posizione dei blocchi nell'isoforma
 * in coordinate della finestra di visualizzazione
 */
function isoform_range(reg) {
	
	var x = d3.scale.log()
				.rangeRound([0, width - 50], .1);
				
	//valorei minimo e massimo di inizio e fine dei blocchi
	min = d3.min(reg, function(d) { return d.start; });
	max = d3.max(reg, function(d) { return d.end; });
	
	x.domain([min, max], .1);
	
	return x;
}

/* SET_SVG
 * 
 * Crea un elemento "svg" dove appendere le gli elementi che 
 * descrivono la struttura dell'isoforma
 */
function set_svg(c, w, h){
	
	
	
	svg = d3.select("body").append("svg")
		.attr("id", c)
		.attr("width", w)
		.attr("height", h);
	
	
	return svg;
}



      
/* DRAW_EXONS
 * box -> variabile che contiente l'elemento "svg"
 * exons -> struttura dati degli esoni
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Disegna la struttura degli esoni. 
 */
function draw_exons(box, exons, x_scale){
	
	console.log(exons);
	  	
	var tip_info = d3.tip()
  			.attr('class', 'd3-tip')
  			.offset([10, 0])
  			.direction("s");
	  			  	
	tf = d3.svg.transform()
		.translate(function (d) { return [x_scale(d.start), 0]; });
		//.scale(0.8, 1);
	
	rect_exons = box.append("g")
		.attr("id", "exons")
		.attr("render-order", 0)
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top + ")");
	//aggiunge i blocchi
	rect_exons.selectAll("rect")
		.data(exons)
		.enter()
		.append("rect")
		.attr("class", function(d) {if(d.alternative == false)
										return "tex";
									else
										return "no_tex"; })
		.attr("width", function(d) { return x_scale(d.end) - x_scale(d.start) - 1; })
		.attr("height", "75")
		.style("stroke", "black")
		.style("fill", function() { return d3.rgb("#228B22"); })
		.style("opacity", 1.0)
		.attr("transform",tf)
		.call(tip_info)
		/*.on("click", function(d){ coord = d3.mouse(this);
								  console.log(coord);
								  if((coord[0] > x_scale(d.start)) | (coord[0] < x_scale(d.end)))
								  		console.log(d.id); })*/
		.on("mouseover", function() { d3.select(this).style('cursor', 'crosshair')
									  .append(tipBlocks.show); })
		.on("mouseout", function() { d3.select(this).style('cursor', 'default')
									 .call(tipBlocks.hide); });
	
			
	return rect_exons;
	
}

/* DRAW_INTRONS
 * box -> variabile che contiente l'elemento "svg"
 * introns -> struttura dati degli introni
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Disegna la struttura degli introni.
 */
function draw_introns(box, introns, x_scale){
	
	line_introns = box.append("g")
		.attr("id", "introns")
		.attr("render-order", -1)
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top + ")");
		
	line_introns.selectAll("line")
		.data(introns)
		.enter().append("line")
		.attr("x1", function(d) { return x(d.start); })
		.attr("y1", 35)
		.attr("x2", function(d) { return x(d.end); })
		.attr("y2", 35)
		.style("stroke", "black")
		.style("stroke-width", 6);
			
	return line_introns;
}



/* CLONE_TRIANGLE_UP
 * svg -> variabile che contiene l'elemento "svg"
 * obj -> oggetto da clonare
 * 
 * Clona l'oggetto contenuto in "obj" e lo aggiunge alla variabile "svg".
 * Il contenitore "g" del segnale alto della tipologia degli splice sites 
 * viene clonato e riutilizzato per aggiungere un segnale basso.
 */
function clone_triangle_up(svg, obj) {
	
	//console.log(obj);
	//console.log(obj.attr("id"));
	
	triangle_down = svg.append("use")
    	.attr("xlink:href","#" + obj.attr("id"))
    	.attr("transform", "translate(0, 140)")
    	.attr("id", "triangle_down");
    
    return triangle_down;
   
}

/* DRAW_SPLICE_SITES
 * box -> variabile che contiente l'elemento "svg"
 * s_s -> struttura dati degli splice sites
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Disegna gli splice sites e i segnali che ne indicano la tipologia.
 */
function draw_splice_sites(box, s_s, x_scale){
	
	tr = d3.svg.symbol()
		.type('triangle-up')
		.size(20);
					
	//aggiunge le splice site
	splice_sites = box.append("g")
		.attr("id", "splice_sites")
		.attr("render-order", 1)
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top + ")");
		
	splice_sites.selectAll("line")
		.data(s_s)
		.enter().append("line")
		.attr("x1", function(d) { if(d.position != null) 
									return x_scale(d.position); })
		.attr("y1", -30)
		.attr("x2", function(d) { if(d.position != null)
									return x_scale(d.position); })
		.attr("y2", 110)
		.style("stroke", "black")
		//.style("stroke-width", 3)
		.style("stroke-dasharray", function(d) { if ((d.type == "init") | (d.type == "term"))
													return 4;
												 else
												 	if(d.type != "unknow")
														return 0; });
	//segnale alto tipologia splice_site		
	triangle_up = box.append("g")
		.attr("id", "triangle_up")
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top + ")");
											
	triangle_up.selectAll("path")
		.data(s_s)
		.enter().append("path")
		.attr("id", function(d, i) { return "up" + i; })
		.attr("d",tr)
		.attr("fill", "red")
		.attr("stroke","#000")
		.attr("stroke-width",1)
		.attr("transform", function (d) {if((d.type == 3) | (d.type == "term"))
											return "translate(" + (x_scale(d.position) - 3) + ",-30)" + "rotate(-90)";
										 else
										 	if((d.type != "unknow") | (d.type == "init"))
												return "translate(" + (x_scale(d.position) + 3) + ",-30)" + "rotate(90)"; });
	
	//viene clonato triangle_up										 
	triangle_down = clone_triangle_up(box, triangle_up);
			
	return splice_sites;
}

/* BUTTON_EXPAND
 * r -> esoni
 * i -> introni
 * s_s -> splice sites
 * tr_d -> segnale della tipologia degli splice sites
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Aggiunge un bottone per l'espansione di tutta la struttura dell'isoforma.
 */
function button_expand(r, i, s_s, tr_d, x_scale){
	
	bt = d3.select("body").append("button")
		.attr("id", "expand_structure")
		.text("Expand")
		.style("top", "20px")
		.style("left", "30px")
		.style("position", "absolute");
			
	bt.on("click", function () {
		
	r.selectAll("rect")
		.transition()
		.duration(750)
		.attr("transform", function(d, i) { return "translate(" + (x_scale(d.start) -1) + "," + i*75 + ")"; });
			
	i.selectAll("line")
		.transition()
		.duration(750)
		.attr("transform", function(d, i) { return "translate(0," + i*45 + ")"; });  
				
	s_s.selectAll("line")
		.transition()
		.duration(750)
		.attr("y2", height);
	
	tr_d.transition()
		.duration(750)
		.attr("transform", function() { return "translate(0," + (height - 100) + ")"; });
		
	});					
}

function button_expand_exons(r, x_scale){
	
	bt = d3.select("body").append("button")
		.attr("id", "expand_structure")
		.text("Expand")
		.style("top", "20px")
		.style("left", "30px")
		.style("position", "absolute");
			
	bt.on("click", function () {
		
	r.selectAll("rect")
		.transition()
		.duration(750)
		.attr("transform", function(d, i) { return "translate(" + (x_scale(d.start) -1) + "," + i*25 + ")"; });
	});
}

/* REGIONS_SCALED
 * r -> regioni
 * 
 * Modifica le regioni in base alla loro dimensione (end - start). Aumenta
 * quelle minori di 50 e diminuisce quelle maggiori di 1500
 */
function regions_scaled(r){
	
	for(i = 0; i < r.length - 1; i++){
		size_regions = r[i].end - r[i].start;
		console.log(size_regions);
		if(size_regions < 50){
			size_regions_scaled = size_regions * (80 / size_regions);
			r[i].end = r[i].end + size_regions_scaled;
			r[i + 1].start = r[i].end;
			console.log(r[i].id);
		}
		if(size_regions > 1500){
			size_regions_scaled = size_regions * (800 / size_regions);
			r[i].end = r[i].end - size_regions_scaled;
			r[i + 1].start = r[i].end;
			console.log(r[i].id);
		}	
	}
			
	return r;
}

//visualizza le informazioni dei blocchi			
tipBlocks = d3.tip()
  			.attr('class', 'd3-tip')
  			.offset([10, 0])
  			.direction("e")
  			.html(function(d) {
   					return "<strong>Block:</strong> <span style='color:red'>" + 
   							d.id + "</span>" + 
   							"<br><strong>Sequence:</strong> <span style='color:green'>" + 
   							d.sequence + "</span>" +
   							"<br><strong>Regions:</strong> <span style='color:yellow'>" + 
   							d.regions + "</span>";
  				});
  				
function pattern_exons(){
	
	defs = d3.select("body").append("svg")
		.append('defs');
	pat_tex = defs.append('pattern')
    	.attr('id', 'pattern_stripe')
    	.attr('patternUnits', 'userSpaceOnUse')
    	.attr('patternTransform', 'rotate(45)')
    	.attr('width', 4)
    	.attr('height', 4)
  		//grandezza stripes
  		.append('rect')
    	.attr('width', 3)
    	.attr('height', 6)
    	.attr('transform', 'translate(0,0)')
    	.style('fill', 'green')
    	.style('opacity', 1.0);
    defs.append('mask')
    	.attr('id', 'mask_stripe')
    	.append('rect')
    	.attr('x', 0)
    	.attr('y', 0)
    	.attr('width', '130%')
    	.attr('height', '130%')
    	.style('fill', 'url(#pattern_stripe)');
    	
    	
   	pat_no_tex = defs.append('pattern')
    	.attr('id', 'pattern_no_stripe')
    	.attr('patternUnits', 'userSpaceOnUse')
    	.attr('patternTransform', 'rotate(45)')
    	.attr('width', 10)
    	.attr('height', 10)
  		.append('rect')
    	.attr('width', '100%')
    	.attr('height', '100%')
    	.attr('transform', 'translate(0,0)')
    	.style('fill', 'green')
    	.style('opacity', 1.0);
    defs.append('mask')
    	.attr('id', 'mask_no_stripe')
    	.append('rect')
    	.attr('x', 0)
    	.attr('y', 0)
    	.attr('width', '130%')
    	.attr('height', '130%')
    	.style('fill', 'url(#pattern_no_stripe)');
}

//carica i dati contenuti nel file json e richiama le funzioni per disegnare la struttura
//dell'isoforma
d3.json("ATP6AP1example2.json", function(error, atp) {
	
	console.log(error);
	isoform = atp[0];
	//console.log(isoform);	
	x = isoform_range(isoform.regions);
	
	regions = regions_scaled(isoform.regions);
	
	
	//console.log(regions);
	boundaries = isoform.boundaries;
	//console.log(boundaries);
	exons = isoform.exons;
	//console.log(exons);
	introns = isoform.introns;
	//console.log(introns);
	
	
	exons_restruct = exons_structure(exons, regions, boundaries);
	introns_restruct = introns_structure(introns, regions, boundaries);
	s_s_restruct = splice_site_structure(boundaries, regions);
	//console.log(exons_restruct);

	//console.log(exons_restruct);
	//console.log(introns_restruct);
	//console.log(s_s_restruct);
	
	svg_box = set_svg("isoform", width + margin_isoform.right + margin_isoform.left, (height - 100) + margin_isoform.top + margin_isoform.bottom);
	svg_box.call(tipBlocks);
	
	pattern_exons();
	
	line_i = draw_introns(svg_box, introns_restruct, x);
	
	s_s = draw_splice_sites(svg_box, s_s_restruct, x);
	
	rect = draw_exons(svg_box, exons_restruct, x);
	
	
	//button_expand(rect, line_i, s_s, triangle_down, x);
	
	button_expand_exons(rect, x);
	//console.log(exons_restruct);
	//console.log(introns_restruct);
	//console.log(s_s_restruct);
	
    
});

	

	