
//margini
var margin_isoform = {top: 100, right: 15, bottom: 15, left: 10};
//dimensione della finestra di visualizzazione dell'isoforma
var height = window.innerHeight + 100 - margin_isoform.top - margin_isoform.bottom,
width = window.innerWidth - margin_isoform.left - margin_isoform.right;


/* EXONS_STRUCTURE
 * extract_exons -> dati degli esoni estratti dal file json
 * extract_regions -> dati delle regioni estratti dal file json
 * extract_boundary -> dati dei boundaries estratti dal file json
 * 
 * Crea la struttura dei dati per gli esoni. Concatena le sequenze
 * nucleotidiche delle regioni. Riporta anche gli ID delle regioni
 * che compongono l'esone.
 */
function exons_structure(extract_exons, extract_regions, extract_boundary){
	
	var l = extract_exons.length;
	
	//array di ogetti "esone"
	var exons = [];
	
	for(i = 0; i < l; i++){
		
		var reg = [];
		
		//left & right boundary
		var l_b = extract_exons[i].left_boundary;
		var r_b = extract_exons[i].right_boundary;
		
		//sequenza nucleotidica
		var seq = "";
		var flag_seq = false;
		var flag_alt = false;
		
		//regione boundary di sinistra
		var region_left = extract_regions[extract_boundary[l_b].first + 1];
		//regione boundary di destra
		var region_right = extract_regions[extract_boundary[r_b].first];
			
		//tipologia regione di sinistra	
		var type_region_left = region_left.type;
		//tipologia regione di destra
		var type_region_right = region_right.type;
		
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
 * extract_boundary -> dati dei boundaries estratti dal file json
 * 
 * Crea la struttura dei dati per gli introni. Concatena le sequenze
 * nucleotidiche delle regioni. Riporta anche gli ID delle regioni
 * che compongono l'introne. Calcola il pattern dai suffissi e prefissi
 * delle regioni che compongono l'introne.
 */
function introns_structure(extract_introns, extract_regions, extract_boundary){
	
	var l = extract_introns.length;
	
	//array di oggetti "introne"
	var introns = [];
	var reg = [];

	
	for(i = 0; i < l; i++){
		
		var reg = [];
		var seq = "";
		var flag_seq = false;
		var flag_intron_ok = true;
		
		//left & right boundary
		var l_b = extract_introns[i].left_boundary;
		var r_b = extract_introns[i].right_boundary;
		
		//regione boundary di sinistra
		var region_left = extract_regions[extract_boundary[l_b].first + 1];
		//regione boundary di destra
		var region_right = extract_regions[extract_boundary[r_b].first];
			
		//tipologia regione di sinistra	
		var type_region_left = region_left.type;
		var alternative_region_left = region_left.alternative;
		//tipologia regione di destra
		var type_region_right = region_right.type;
		var alternative_region_right = region_right.alternative;
		
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
		
		    //suffisso e prefisso della sequenza nucleotidica
			var l_suffix = extract_introns[i].suffix.length;
			var pattern = extract_introns[i].prefix.substr(0, 2).concat(extract_introns[i].suffix.substr(l_suffix - 2, l_suffix));
			
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
	
	var l = extract_boundaries.length;
	
	//array di oggetti "splice sites"
	var s_s = [];
	
	for(i = 0; i < l; i++){
		if(extract_boundaries[i].first == -1){
			pos_b = null;
			t = "unknow";
		}	
		else{
			t = extract_boundaries[i].type;	
			if((t == "5") | (t == "both"))
				pos_b = extract_regions[extract_boundaries[i].first + 1].start;
			if((t == "3") | (t == "term"))
				pos_b = extract_regions[extract_boundaries[i].first].end;
			if(t == "init")
				pos_b = extract_regions[extract_boundaries[i].first].start;
						
			s_s.push({
				"position" : pos_b,
				"type" : t			
			});
		}
	}	
	return s_s;
}

/* ISOFORM_SCALE
 * reg -> dati delle regioni estratti dal file json
 * 
 * Definisce un range per trasformare la posizione dei blocchi nell'isoforma
 * in coordinate della finestra di visualizzazione
 */
function isoform_range(reg) {
	
	var x = d3.scale.log()
				.rangeRound([0, width - 50], .1);
				
	//valorei minimo e massimo di inizio e fine dei blocchi
	var min = d3.min(reg, function(d) { return d.start; });
	var max = d3.max(reg, function(d) { return d.end; });
	
	x.domain([min, max], .1);
	
	return x;
}

/* SET_SVG
 * c -> classe
 * w -> width
 * h -> height
 * p -> array per la posizione
 * 
 * Crea un elemento "svg", specificando la posizione, 
 * la dimensione e la classe
 */
function set_svg(c, w, h, p){
	
	var svg = d3.select("body").append("svg")
		.attr("id", c)
		.attr("width", w)
		.attr("height", h)
		.style("position", p[0])
		.style("left", p[1])
		.style("right", p[2])
		.style("top", p[3]);	
	
	return svg;
}



//----------------FUNZIONI PER GLI ESONI---------------------------------------------------------

//----------------FUNZIONI PER LA FINESTRA INFO--------------------------------------------------


/* WINDOWS_INFO_SCALE
 * reg -> regioni dell'esone selezionato
 * h_info -> altezza finestra 
 * 
 * Ritorna una funzione che riscala ogni valore nel range 
 * della finestra di visualizzazione
 */
function window_info_scale(reg, h_info){
    
    var y = d3.scale.log()
              .rangeRound([0, h_info/2], .1);
                
    //valorei minimo e massimo di inizio e fine dei blocchi
    var min_r = d3.min(reg[0], function(d) { return d.start; });
    var max_r = d3.max(reg[0], function(d) { return d.end; });
    
    if(reg.length > 1){
        min_i = d3.min(reg[1], function(d) { return d.start; });
        max_i = d3.max(reg[1], function(d) { return d.end; });
    }
    else{
        min_i = min_r;
        max_i = max_r;
    }
    
    if(min_r <= min_i)
        min = min_r;
    else
        min = min_i; 
    
    if(max_r >= max_i)
        max = max_r;
    else
        max = max_i;
    
    y.domain([min, max], .1);
    
    return y;
}

/* SVG_INFO_BOX
 * 
 * Crea una finestra dove saranno visualizzati gli elementi
 * appartenenti alla selezione. Viene aggiunto il tasto per
 * fare un clean della finestra e riattivare la struttura del gene
 */
function svg_info_box(){
    
    s_w = 650;
    s_h = 300;
                
    var p_s = ["absolute", "20px", "10px", "480px", "10px"];
                                   
    var s_i = set_svg("expande_info", s_w, s_h, p_s);
    
    s_i.attr("viewbox", "0 0 500 400")
        .style("border", "3px solid #cccccc");
        
    d3.select("body").append("button")
        .attr("id", "clear_vis")
        .text("Clear")
        .style("top", "456px")
        .style("left", "619px")
        .style("position", "absolute")
        .on("click", function(){ 
                        d3.select("#regions_selected").remove(); 
                        
                        d3.selectAll("#exon")
                           .attr("pointer-events", "yes")
                           .style("stroke-width", 0)
                           .style("fill", function() { return d3.rgb("#228B22"); });
                                   
                        d3.selectAll("#exon_stripe")
                          .attr("pointer-events", "yes")
                          .style("stroke-width", 0)
                          .style("fill", 'url(#diagonalHatch)');
                        
                        d3.selectAll("#exon")
                          .attr("pointer-events", "yes")
                          .style("stroke-width", 0)
                          .style("fill", function() { return d3.rgb("#228B22"); });
                        
                        d3.selectAll("#intron")
                            .style("stroke", "black")
                            .attr("pointer-events", "yes"); }); 
    return s_i;
}

/* DISPLAY_INFO
 * s_i -> finestra creata per visualizzare le informazioni
 * domain -> dominio della finestra di visualizzazione degli elementi
 * elements -> elementi estratti dalla struttura del gene
 *             in base alla selezione.
 * 
 * Visualizza gli elementi appartenenti alla selezione.
 */
function display_info(s_i, domain, elements){
    
    var exons_info = elements[0];
    var introns_info = elements[1];
    
    var tipElement = d3.tip()
        .attr('class', 'd3-tip')
        .offset([10, 0])
        .direction("e")
        .html(function(d) {
                return "<strong>Region:</strong> <span style='color:red'>" + 
                        d.id + "</span>" + 
                        "<br><strong>Start:</strong> <span style='color:yellow'>" + 
                        d.start + "</span>" +
                        "<br><strong>End:</strong> <span style='color:yellow'>" + 
                        d.end + "</span>" +
                        "<br><strong>Sequence:</strong> <span style='color:green'>" + 
                        d.sequence + "</span>";
                    
        });
        
    var tf_info_ex = d3.svg.transform()
        .translate(function (d, i) { return [15, i * 50]; })
        .scale(function (d, i) { return [2.8, 1]; });
    
    var tf_info_in = d3.svg.transform()
        .translate(function (d, i) { return [15, (i * 50 + exons_info.length * 50)]; })
        .scale(function (d, i) { return [2, 1]; });
        
    
    var tf_g = d3.svg.transform()
        .translate(function (d, i) { return [15, 20]; });
    
    var g = s_i.append("g")
        .attr("id", "regions_selected")
        .attr("transform", tf_g)
        .call(tipElement);
        
    g.selectAll("rect")
        .data(exons_info)
        .enter().append("rect")
        .attr("width", function(d) { return domain(d.end) - domain(d.start); })
        .attr("height", 40)
        .style("fill", function() { return d3.rgb("#B8860B"); })
        .style("opacity", 0.7)
        .attr("transform", tf_info_ex)
        .on("mouseover", tipElement.show)
        .on("mouseout", tipElement.hide);
    
    if(introns_info != null){
    	g.selectAll("line")
    		.data(introns_info)
			.enter().append("line")
			.attr("x1", function(d) { return domain(d.start); })
			.attr("y1", 35)
			.attr("x2", function(d) { return domain(d.end); })
			.attr("y2", 35)
			.attr("transform", tf_info_in)
			.style("stroke", "black")
			.style("stroke-width", 6);   
	}     
}

/* DISPLAY_INFO_STRIPE
 * s_i -> finestra creata per visualizzare le informazioni
 * domain -> dominio della finestra di visualizzazione degli elementi
 * elements -> elementi estratti dalla struttura del gene
 *             in base alla selezione.
 * 
 * Visualizza gli elementi appartenenti alla selezione 
 * (per esone 'alternative').
 */     
function display_info_stripe(s_i, domain, elements){
    
    var exons_info = elements[0];
    var introns_info = elements[1];
    
    var tipElement = d3.tip()
        .attr('class', 'd3-tip')
        .offset([10, 0])
        .direction("e")
        .html(function(d) {
                return "<strong>Region:</strong> <span style='color:red'>" + 
                        d.id + "</span>" + 
                        "<br><strong>Start:</strong> <span style='color:yellow'>" + 
                        d.start + "</span>" +
                        "<br><strong>End:</strong> <span style='color:yellow'>" + 
                        d.end + "</span>" +
                        "<br><strong>Sequence:</strong> <span style='color:green'>" + 
                        d.sequence + "</span>";
                    
        });
        
    var tf_info_ex = d3.svg.transform()
        .translate(function (d, i) { return [15, i * 50]; })
        .scale(function (d, i) { return [2.8, 1]; });
    
    var tf_info_in = d3.svg.transform()
        .translate(function (d, i) { return [15, (i * 50 + exons_info.length * 50)]; })
        .scale(function (d, i) { return [2, 1]; });
    
    var tf_g = d3.svg.transform()
        .translate(function (d, i) { return [15, 20]; });
    
    var g = s_i.append("g")
        .attr("id", "regions_selected")
        .attr("transform", tf_g)
        .call(tipElement);
        
    g.selectAll("rect")
        .data(exons_info)
        .enter().append("rect")
        .attr("width", function(d) { return domain(d.end) - domain(d.start); })
        .attr("height", 40)
        .style("fill", function() { return d3.rgb("#B8860B"); })
        .style("opacity", 0.7)
        .attr("transform", tf_info_ex)
        .on("mouseover", tipElement.show)
        .on("mouseout", tipElement.hide);
        
    if(introns_info != null){
    	g.selectAll("line")
    		.data(introns_info)
			.enter().append("line")
			.attr("x1", function(d) { return domain(d.start); })
			.attr("y1", 35)
			.attr("x2", function(d) { return domain(d.end); })
			.attr("y2", 35)
			.attr("transform", tf_info_in)
			.style("stroke", "black")
			.style("stroke-width", 6);   
	}     
}

//variabile che contiente la funzione per 
//visualizzare le informazioni sugli esoni  
//DA SISTEMARE        
tipBlocks = d3.tip()
    .attr('class', 'd3-tip')
    .offset([10, 0])
    .direction("e")
    .html(function(d) {
            return "<strong>Esone:</strong> <span style='color:red'>" + 
                    d.id + "</span>" + 
                    "<br><strong>Sequence:</strong> <span style='color:green'>" + 
                    d.sequence + "</span>" +
                    "<br><strong>Regions:</strong> <span style='color:yellow'>" + 
                    d.regions + "</span>";
        });


/* PATTERN_EXONS
 * 
 * Definisce un pattern a strisce per differenziare la tipologia
 * degli esoni.
 */                
function pattern_exons(){
    
    defs = d3.select("body").append("svg")
        .append('defs');
    defs.append('defs')
        .append('pattern')
        .attr('id', 'diagonalHatch')
        .attr('patternUnits', 'userSpaceOnUse')
        .attr('width', 4)
        .attr('height', 4)
        .append('path')
        .attr('d', 'M-1,1 l2,-2 M0,4 l4,-4 M3,5 l2,-2')
        .attr('stroke', '#000000')
        .attr('stroke-width', 1.5);
}


/* CHECK_STRUCTURE_ELEMENT
 * regions_info -> vettore che contiene gli ID delle regioni
 *                 che compongono l'esone 
 * c -> coordinate(mouse) della selezione
 * 
 * Crea una struttura di esoni(per regioni) e introni in base
 * alle coordinate della selezione.
 */
function check_structure_element(regions_info, c){
                          
    var r_info = [];
    var i_info = [];
    var c_info = [];
    
    for(rg = 0; rg < regions_info.length; rg++)
        r_info.push(regions[regions_info[rg]]);
    
    for(i = 0; i < introns_restruct.length; i++)
        if((c >= introns_restruct[i].start) & (c <= introns_restruct[i].end))
            i_info.push(introns_restruct[i]);
            
    if(r_info.length != 0)
        c_info.push(r_info);
    if(i_info.length != 0)
        c_info.push(i_info);
    
    return c_info;   
}

/* DRAW_EXONS
 * box -> variabile che contiente l'elemento "svg"
 * exons -> struttura dati degli esoni
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Disegna gli esoni differenziandoli in base al campo "alternative". Aggiunge
 * i pop-up per le informazioni e chiama le funzioni per la selezione degli 
 * elementi.
 */
function draw_exons(box, exons, x_scale){
	
	var exons_stripe = [];
	for(k = 0; k < exons.length; k++)
	   if(exons[k].alternative == false)
	       exons_stripe.push(exons[k]);
    	  			  	
	var tf = d3.svg.transform()
		.translate(function (d) { return [x_scale(d.start), 0]; });
		//.scale(0.8, 1);
	
	var rect_exons = box.append("g")
		.attr("id", "exons")
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top + ")");
	//aggiunge i blocchi
	rect_exons.selectAll("rect")
		.data(exons)
		.enter()
		.append("rect")
		.attr("id", "exon")
		.attr("width", function(d) { return x_scale(d.end) - x_scale(d.start) - 1; })
		.attr("height", "0")
		//.style("stroke", "black")
		.style("fill", function() { return d3.rgb("#228B22"); })
		.style("opacity", 1.0)
		.attr("transform",tf)
		.on("click", function(d){ 
		               if(d.alternative == true){
		                  d3.selectAll("#exon")
						    .style("fill", function() { return d3.rgb("#808080"); })
							.attr("pointer-events", "none");
							
						  d3.selectAll("#intron")
                            .style("stroke", function() { return d3.rgb("#808080"); })
                            .attr("pointer-events", "none");
									
						  d3.select(this)
							.style("fill", function() { return d3.rgb("#228B22"); })
							.style("stroke", function() { return d3.rgb("#B8860B"); })
							.style("stroke-width", 3);
						  
						  var coord_x = x_scale.invert(d3.event.pageX);
						  
						  if((coord_x > x_scale(d.start)) | (coord_x < x_scale(d.end))){
						      info_structure = check_structure_element(d.regions, coord_x);  	
						      x_info = window_info_scale(info_structure, s_h);
						      display_info(svg_info, x_info, info_structure);}
        			      }})							  		
		.on("mouseover", function() { d3.select(this).style('cursor', 'crosshair')
									  .append(tipBlocks.show); })
		.on("mouseout", function() { d3.select(this).style('cursor', 'default')
									 .call(tipBlocks.hide); })
	    .transition()
	    .duration(750)
	    .attr("height", 75);
	
	//esoni alternative = false								 
	var rect_exons_stripe = box.append("g")
        .attr("id", "exons")
        .attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top + ")");
    //aggiunge i blocchi
    rect_exons_stripe.selectAll("rect")
        .data(exons_stripe)
        .enter()
        .append("rect")
        .attr("id", "exon_stripe")
        .attr("width", function(d) { return x_scale(d.end) - x_scale(d.start) - 1; })
        .attr("height", "0")
        //.style("stroke", "black")
        .style("fill", 'url(#diagonalHatch)')
        .attr("transform",tf) 
        .on("click", function(d){ 
                        d3.selectAll("#exon_stripe")
                          .style("fill", function() { return d3.rgb("#808080"); })
                          .attr("pointer-events", "none");
                        d3.selectAll("#exon")
                          .style("fill", function() { return d3.rgb("#808080"); })
                          .attr("pointer-events", "none");
                        d3.selectAll("#intron")
                            .style("stroke", function() { return d3.rgb("#808080"); })
                            .attr("pointer-events", "none");
                                    
                        d3.select(this)
                          .style("fill", function() { return d3.rgb("#228B22"); })
                          .style("stroke", function() { return d3.rgb("#B8860B"); })
                          .style("stroke-width", 3);
                        var coord_x = x_scale.invert(d3.event.pageX);
                          //console.log(coord);
                        if((coord_x > x_scale(d.start)) | (coord_x < x_scale(d.end))){
                              //console.log(d.id); 
                              //console.log(d.regions);    
                            info_structure = check_structure_element(d.regions, coord_x);  
                              
                            x_info = window_info_scale(info_structure, s_h);
                            display_info_stripe(svg_info, x_info, info_structure);}})                           
        .on("mouseover", function() { d3.select(this).style('cursor', 'crosshair')
                                      .append(tipBlocks.show); })
        .on("mouseout", function() { d3.select(this).style('cursor', 'default')
                                     .call(tipBlocks.hide); })
        .transition()
        .duration(750)
        .attr("height", 75);
        		
	return rect_exons;	
}

//----------------------------------------------------------------------------------------------



//------------------------------------FUNZIONI PER GLI INTRONI----------------------------------

/* DRAW_INTRONS
 * box -> variabile che contiente l'elemento "svg"
 * introns -> struttura dati degli introni
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Disegna gli introni.
 */
function draw_introns(box, introns, x_scale){
    	
	var line_introns = box.append("g")
		.attr("id", "introns")
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top + ")");
		
	line_introns.selectAll("line")
		.data(introns)
		.enter().append("line")
		.attr("id", "intron")
		.attr("x1", function(d) { return x(d.start); })
		.attr("y1", 35)
		.attr("x2", function(d) { return x(d.start); })
		.attr("y2", 35)
		.style("stroke", "black")
		.style("stroke-width", 6)
		.transition()
		.delay(750)
		.duration(750)
		.attr("x2", function(d) { return x(d.end); });
			
	return line_introns;
}

//----------------------------------------------------------------------------------------------

//-------------------------------------FUNZIONI PER GLI SPLICE_SITES----------------------------

/* CLONE_TRIANGLE_UP
 * svg -> variabile che contiene l'elemento "svg"
 * obj -> oggetto da clonare
 * 
 * Clona l'oggetto contenuto in "obj" e lo aggiunge alla finestra di 
 * visualizzazione contenuta nella variabile "svg".
 * Il contenitore "g" del "segnale alto" della tipologia degli splice sites 
 * viene clonato e riutilizzato per aggiungere un "segnale basso".
 */
function clone_triangle_up(svg, obj) {
	
	var triangle_down = svg.append("use")
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
	
	//variabile per i simboli della tipologia 
	//di splice_sites
	var tr = d3.svg.symbol()
		.type('triangle-up')
		.size(20);
					
	//aggiunge le splice site
	var splice_sites = box.append("g")
		.attr("id", "splice_sites")
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
		.style("opacity", "0.0")
		.style("stroke", "black")
		//.style("stroke-width", 3)
		.style("stroke-dasharray", function(d) { 
		                              if ((d.type == "init") | (d.type == "term"))
									       return 4;
									  else
										   if(d.type != "unknow")
										      return 0; })
	    .transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
										  
	//segnale alto tipologia splice_site		
	var triangle_up = box.append("g")
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
		.style("opacity", "0.0")
		.attr("transform", function (d) {
		                      if((d.type == 3) | (d.type == "term"))
							     return "translate(" + (x_scale(d.position) - 3) + ",-30)" + "rotate(-90)";
							  else
							     if((d.type != "unknow") | (d.type == "init"))
								    return "translate(" + (x_scale(d.position) + 3) + ",-30)" + "rotate(90)"; })
	    .transition()
        .delay(1500)
        .duration(750)
        .style("opacity", "1.0");
        
	//viene clonato triangle_up										 
	var triangle_down = clone_triangle_up(box, triangle_up);
			
	return splice_sites;
}

//--------------------------------------------------------------------------------------------------------------


/* REGIONS_SCALED
 * r -> regioni
 * 
 * Modifica le regioni in base alla loro dimensione (size = end - start). Aumenta
 * quelle minori di 50 e diminuisce quelle maggiori di 1500
 */
function regions_scaled(r){
	
	for(i = 0; i < r.length - 1; i++){
		size_regions = r[i].end - r[i].start;
		
		//regioni troppo piccole (< 50 pixel)
		if(size_regions < 50){
			size_regions_scaled = size_regions * (80 / size_regions);
			r[i].end = r[i].end + size_regions_scaled;
			r[i + 1].start = r[i].end;
		}
		//regioni troppo grandi (> 1500 pixel)
		if(size_regions > 1500){
			size_regions_scaled = size_regions * (800 / size_regions);
			r[i].end = r[i].end - size_regions_scaled;
			r[i + 1].start = r[i].end;
		}	
	}
    return r;
}

/* CHANGE_GENE
 * 
 * Permette di cambiare gene, rimuovendo la struttura disegnata e 
 * richiamando nuovamente la funzione "init" per disegnare la nuova
 * struttura.
 */
function change_gene(){
    
    //rimozione della struttura
    var g = d3.select("#isoform").selectAll("g");
    g.remove();
    
    //rimozione del titolo
    d3.select("#title").remove();
    
    init();   
}

/* SELECT_GENE
 * 
 * Crea il menù per selezionare il gene di cui visualizzare la struttura.
 * Inizializza le finestre di visualizzazione e chiama la funzione "init"
 */
function select_gene(){
      
    var s_g = d3.select("body").append("g");
    s_g.style("top", "15px")
       .style("left", "25px")
       .style("position", "absolute");
    s_g.html(function() { return '<select id="selezione"><option value="ATP6AP1example2">ATP6AP1</option>' + 
                                   '<option value="ATP6AP1example3">ATP6AP1_v2</option>' + 
                                   '</select>'; }); 
    
    s = d3.select("#selezione")
            .on("change", change_gene);
    
    s.property("value", "ATP6AP1example2");
    
    
    //texture esoni conservativi
    pattern_exons();
    
    pos_box = ["absolute", "20px", "10px", "110px", "10px"];
    svg_box = set_svg("isoform", width - margin_isoform.right + margin_isoform.left, 300, pos_box);
    svg_box.style("border", "3px solid #cccccc")
         .call(tipBlocks);
        
    init();   
    
    svg_info = svg_info_box();                           
}


/* DISPLAY_GENE
 * id -> stringa che identifica il gene da visualizzare
 * 
 * Visualizza il nome del gene di cui sarà visualizzata
 * la struttura.
 */
function display_gene(id){
    
    var title = "";
    
    title = title.concat(id + " gene structure");
    
    var pos_title = ["absolute", "20px", "10px", "40px", "10px"];
    
    var svg_title = set_svg("title", width - margin_isoform.right + margin_isoform.left, 50, pos_title);
    svg_title.style("border", "1px solid #cccccc");
    svg_title.append("text")
       .attr("id", "title")
       .attr("x", 10)
       .attr("y", 35)
       .attr("font-family", "Arial, Helvetica, sans-serif")
       .attr("font-size", "30px")
       .style("fill", "black")
       .style("opacity", "0.0")
       .text(title)
       .transition()
       .duration(1000)
       .style("opacity", "1.0");
}

/* COPY_STRUCTURE
 * s -> struttura da copiare
 * 
 * Copia il vettore delle regioni che poi sara modificato
 * per scalare le dimensioni degli elementi
 */
function copy_structure(s){
    
    var copy_reg = [];
    for(t = 0; t < s.length; t++)
        copy_reg.push({
            "start" : s[t].start,
            "end" : s[t].end,
            "sequence" : s[t].sequence,
            "type" : s[t].type,
            "alternative" : s[t].alternative,
            "coverage" : s[t].coverage,
            "last" : s[t].last,
            "id" : s[t].id    
        });
    
    return copy_reg;    
}

/* INIT
 * 
 * Inizializza tutte le funzione per disegnare la struttura. 
 * Carica i dati dal file json relativo al gene selezionato.
 * Di default carica "ATP6AP1.json"
 */
function init(){
    
    var string = "";
    string = string.concat(s.property("value"), ".json");
    
    //console.log(string);
    //carica i dati contenuti nel file json e richiama le funzioni per disegnare la struttura
    //dell'isoforma
    d3.json(string, function(error, atp) {
	
	   console.log(error);
	   isoform = atp[0];
	   //console.log(isoform);	
	
	   original_regions = copy_structure(isoform.regions);
	   //regioni
	   x = isoform_range(isoform.regions);
	   regions = regions_scaled(isoform.regions);
	
	   display_gene(isoform.gene);
	
	   //boundaries
	   boundaries = isoform.boundaries;
	   //console.log(boundaries);
	   //esoni
	   exons = isoform.exons;
	   //console.log(exons);
	   //introni
	   introns = isoform.introns;
	   //console.log(introns);
	
	   //esoni, introni e boundaries ricostruiti
	   exons_restruct = exons_structure(exons, regions, boundaries);
	   introns_restruct = introns_structure(introns, regions, boundaries);
	   s_s_restruct = splice_site_structure(boundaries, regions);
	   //console.log(exons_restruct);

	   //console.log(exons_restruct);
	   //console.log(introns_restruct);
	   //console.log(s_s_restruct);
	
	   //finestra che contiene la struttura del gene
	   
	
	   line_i = draw_introns(svg_box, introns_restruct, x);
	   rect = draw_exons(svg_box, exons_restruct, x);
	   s_s = draw_splice_sites(svg_box, s_s_restruct, x);
	  	
    });
}







	

	