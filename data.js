
//margini
var margin_isoform = {top: 100, right: 15, bottom: 15, left: 10};
//dimensione della finestra di visualizzazione dell'isoforma
var height = window.innerHeight + 100 - margin_isoform.top - margin_isoform.bottom;
var width = window.innerWidth;
// - margin_isoform.left - margin_isoform.right;


/* EXONS_STRUCTURE
 * extract_exons -> dati degli esoni estratti dal file json
 * extract_regions -> dati delle regioni estratti dal file json
 * extract_boundary -> dati dei boundaries estratti dal file json
 * 
 * Crea la struttura dei dati per gli esoni. Concatena le sequenze
 * nucleotidiche delle regioni. Riporta anche gli ID delle regioni
 * che compongono l'esone.
 */
function exons_structure (extract_exons, extract_regions, extract_boundary) {
	
		
	//array di oggetti "esone"
	var exons = [];
	
	//numero di esoni
	var l = extract_exons.length;
	
	for (i = 0; i < l; i++) {
		
		var reg = [];
		
		//proprietà esone
		var exon_prop = {
			//left & right boundary
			l_b : extract_exons[i].left_boundary,
			r_b : extract_exons[i].right_boundary,
			annot : extract_exons[i].annotated,
			//sequenza nucleotidica
			seq : "",
			flag_seq : false,
			flag_alt : false,
		};
		
		var region_prop = {
			//regione boundary di sinistra
			r_l : extract_regions[extract_boundary[exon_prop.l_b].first + 1],
			//regione boundary di destra
			r_r : extract_regions[extract_boundary[exon_prop.r_b].first]
		};
		
		var boundary_prop = {	
			//tipologia regione di sinistra	
			t_r_l : region_prop.r_l.type,
			//tipologia regione di destra
			t_r_r : region_prop.r_r.type
		};
		
		//esclude le regioni non codificanti
		if ((boundary_prop.t_r_l == "codifying") & (boundary_prop.t_r_l != "unknow")) {
			start_exon = region_prop.r_l.start;
			if ((boundary_prop.t_r_r == "codifying") & (boundary_prop.t_r_r != "unknow"))
				end_exon = region_prop.r_r.end;
			
			//assembla la sequenza nucleotidica e salva le regioni appartenenti all'esone
			for (j = region_prop.r_l.id; j <= region_prop.r_r.id; j++) {
			
				if (extract_regions[j].sequence == null)
					exon_prop_flag.seq = true;
				else {
					if(extract_regions[j].type == "codifying")
						exon_prop.seq = exon_prop.seq.concat(extract_regions[j].sequence);
					else 
						exon_prop.flag_seq = true;
					}
			
				reg.push(extract_regions[j].id);
				//se tutte le regioni sono "alternative" allora
				//anche l'esone è "alternative" CHIEDERE!!
				if(extract_regions[j].alternative == true)
					exon_prop.flag_alt = true;
			}
			
			if (exon_prop.flag_seq == true)
				exon_prop.seq = null;
				
			//costruisce l'oggetto esone
			exons.push({
					"id" : i,
					"start" : start_exon,
					"end" : end_exon,
					"sequence" : exon_prop.seq,
					"regions" : reg,
					"alternative" : exon_prop.flag_alt,
					"annotated" : exon_prop.annot
			});
			reg = [];
		}
	}
	console.log(exons);
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
	
	//numero di introni
	var l = extract_introns.length;
	
	//array di oggetti "introne"
	var introns = [];
	
	for(i = 0; i < l; i++){
		
		var reg = [];
		
		var intron_prop = {
			//stringa per la sequenza nucleotidica
			seq : "",
			flag_seq : false,
			flag_intron_ok : true,
		
			//left & right boundary
			l_b : extract_introns[i].left_boundary,
			r_b : extract_introns[i].right_boundary,
			pattern : ""
		};
		
		var region_prop = {
			//regione boundary di sinistra
			r_l : extract_regions[extract_boundary[intron_prop.l_b].first + 1],
			//regione boundary di destra
			r_r : extract_regions[extract_boundary[intron_prop.r_b].first]
		};
		
		var boundary_prop = {
			//tipologia regione di sinistra	
			t_r_l : region_prop.r_l.type,
			//alternative
			a_r_l : region_prop.r_l.alternative,
			//tipologia regione di destra
			t_r_r : region_prop.r_r.type,
			//alternative
			a_r_r : region_prop.r_r.alternative
		};
		
		//esclude le regioni non codificanti
		if((boundary_prop.t_r_l != "unknow") & (boundary_prop.t_r_r != "unknow")){
			if((extract_boundary[intron_prop.l_b].type ==  "5") | (extract_boundary[intron_prop.l_b].type == "both")){
				if((boundary_prop.t_r_l == "codifying") & (boundary_prop.a_r_l == true))
					start_intron = region_prop.r_l.start;
				else
					if(boundary_prop.t_r_l == "spliced")
						start_intron = region_prop.r_l.start;
					else
						intron_prop.flag_intron_ok == false;
			}
			else
				intron_prop.flag_intron_ok = false;
			
			if((extract_boundary[intron_prop.r_b].type ==  "3") | (extract_boundary[intron_prop.r_b].type == "both")){
				if((boundary_prop.t_r_r == "codifying") & (boundary_prop.a_r_r == true))
					end_intron = region_prop.r_r.start;
				else
					if(boundary_prop.t_r_r == "spliced")
						end_intron = region_prop.r_r.end;
					else
						intron_prop.flag_intron_ok = false;
			}
			else
				intron_prop.flag_intron_ok = false;
		}
		else
			intron_prop.flag_intron_ok = false;
			
		if(intron_prop.flag_intron_ok){
			//assembla la sequenza nucleotidica e salva le regioni appartenenti all'introne
			for(j = region_prop.r_l.id; j <= region_prop.r_r.id; j++){
			
				if(extract_regions[j].sequence == null)
					intron_prop.flag_seq = true;
				else{
					if(((extract_regions[j].type == "codifying") & (extract_regions[j].alternative == true)) | extract_regions[j].type == "spliced")
						intron_prop.seq = intron_prop.seq.concat(extract_regions[j].sequence);
					else
						intron_prop.flag_seq = true;
					}
			
				reg.push(extract_regions[j].id);
			}
			
		    if(intron_prop.flag_seq == true)
				intron_prop.seq = null;
		
		    //suffisso e prefisso della sequenza nucleotidica
			var l_suffix = extract_introns[i].suffix.length;
			intron_prop.pattern = extract_introns[i].prefix.substr(0, 2).concat(extract_introns[i].suffix.substr(l_suffix - 2, l_suffix));
			
			introns.push({
					"start" : start_intron,
					"end" : end_intron,
					"sequence" : intron_prop.seq,
					"suffix" : extract_introns[i].suffix,
					"prefix" : extract_introns[i].prefix,
					"pattern" : intron_prop.pattern,
					"regions" : reg,
					"id" : i
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
	
	//numero dei boundaries
	var l = extract_boundaries.length;
	
	//array di oggetti "splice sites"
	var s_s = [];
	
	var boundary_prop = {
		//posizione
		pos : null,
		//tipologia
		t : ""
	};
	
	for(i = 0; i < l; i++){
		if(extract_boundaries[i].first == -1){
			boundary_prop.pos = null;
			boundary_prop.t = "unknow";
		}	
		else{
			boundary_prop.t = extract_boundaries[i].type;	
			if((boundary_prop.t == "5") | (boundary_prop.t == "both"))
				boundary_prop.pos = extract_regions[extract_boundaries[i].first + 1].start;
			if((boundary_prop.t == "3") | (boundary_prop.t == "term"))
				boundary_prop.pos = extract_regions[extract_boundaries[i].first].end;
			if(boundary_prop.t == "init")
				boundary_prop.pos = extract_regions[extract_boundaries[i].first].start;
						
			s_s.push({
				"position" : boundary_prop.pos,
				"type" : boundary_prop.t			
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
	
	//range per la finestra della struttura
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
		.style("position", p.pos)
		.style("left", p.left)
		.style("right", p.right)
		.style("top", p.top);	
	
	return svg;
}


/* WINDOWS_INFO_SCALE
 * reg -> regioni dell'esone selezionato
 * h_info -> altezza finestra 
 * 
 * Ritorna una funzione che riscala ogni valore nel range 
 * della finestra di visualizzazione
 */
function window_info_scale(reg, h_info){
    
    //range per la finestra degli elementi selezionati
    y = d3.scale.log()
              .rangeRound([0, h_info / 2], .1);
                
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


/* LEGEND_BOX
 * 
 * Visualizza la legenda per gli elementi delle struttura
 * del gene
 */
function legend_box(){

    //dimensioni fisse degli elementi della legenda
    var l_w = window.innerWidth - width / 2 - margin_isoform.left - margin_isoform.right - 5;
    var l_h = 50;
    var off_set = 50;
    var off_set_y = 25;
    //distanziatore delle colonne della legenda
    //deve descrescere proporzionalmente al numero
    //delle colonne
    var c = 40;
    var height_exon = 18;
    
    //vettore per il posizionamento            
    var p_s = {
    	pos : "absolute", 
    	left : (window.innerWidth - l_w - margin_isoform.right) +  "px",
    	right : margin_isoform.right + "px",
    	top: "40px", 
    	bottom : "10px"
    };
    
   
    
    //variabile per traslare gli elementi
	var tf_element = d3.svg.transform()
		.translate(function () { return [70, 3]; });
	var tf_text = d3.svg.transform()
		.translate(function () { return [70, 7]; });
	
	var color_exon = function() { return d3.rgb("#228B22"); };
	var color_exon_stripe = function() { return 'url(#diagonalHatch)'; };
    var color_intron = function() { return d3.rgb("black"); };
    var color_splice_site = function() { return d3.rgb("black"); };
                                        
    var s_l = set_svg("legend", l_w, l_h, p_s);
   
    s_l.attr("viewbox", function() { return "0 0" + l_w + l_h; });
 
    //esoni
    s_l.append("rect")
    	.attr("x", 0)
    	.attr("y", 0)
    	.attr("width", 40)
    	.attr("height", height_exon)
    	.attr("transform", tf_element)
    	.style("fill", color_exon)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    s_l.append("text")
    	.attr("x", off_set)
    	.attr("y", height_exon/2)
    	.style("font-size", "12px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.style("opacity", "0.0")
    	.attr("transform", tf_text)
    	.text("Alternative exons")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    
    s_l.append("rect")
    	.attr("x", 0)
    	.attr("y", off_set_y)
    	.attr("width", 40)
    	.attr("height", height_exon)
    	.attr("transform", tf_element)
    	.style("fill", color_exon)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    s_l.append("rect")
    	.attr("x", 0)
    	.attr("y", off_set_y)
    	.attr("width", 40)
    	.attr("height", height_exon)
    	.attr("transform", tf_element)
    	.style("fill", color_exon_stripe)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    s_l.append("text")
    	.attr("x", off_set)
    	.attr("y", height_exon/2 + off_set_y)
    	.style("font-size", "12px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.style("opacity", "0.0")
    	.attr("transform", tf_text)
    	.text("Conservative exons")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    
    //introni	
    s_l.append("line")
    	.attr("x1", off_set * 4 + c)
    	.attr("y1", 10)
		.attr("x2", off_set * 4 + c + 40)
		.attr("y2", 10)
    	.attr("transform", tf_element)
    	.style("stroke", color_intron)
    	.style("stroke-width", 6)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    s_l.append("text")
    	.attr("x", off_set * 5 + c)
    	.attr("y", 10)
    	.style("font-size", "12px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.style("opacity", "0.0")
    	.attr("transform", tf_text)
    	.text("Introns")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    
    //splice sites
    s_l.append("line")
    	.attr("x1", off_set * 4 + c)
    	.attr("y1", off_set_y)
		.attr("x2", off_set * 4 + c)
		.attr("y2", off_set_y + 20)
    	.attr("transform", tf_element)
    	.style("stroke", color_intron)
    	.style("stroke-width", 2)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    s_l.append("line")
    	.attr("x1", off_set * 4 + c + 10)
    	.attr("y1", off_set_y)
		.attr("x2", off_set * 4 + c + 10)
		.attr("y2", off_set_y + 20)
    	.attr("transform", tf_element)
    	.style("stroke", color_intron)
    	.style("stroke-width", 2)
    	.style("stroke-dasharray", 4)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
    s_l.append("text")
    	.attr("x", off_set * 5 + c)
    	.attr("y", off_set_y + 10)
    	.style("font-size", "12px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.style("opacity", "0.0")
    	.attr("transform", tf_text)
    	.text("Splice sites").transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0"); 
    
    var s_sy = d3.svg.symbol()
		.type('triangle-up')
		.size(10);
	//type 3
    s_l.append("line")
    	.attr("x1", off_set * 8 + c)
    	.attr("y1", 4)
		.attr("x2", off_set * 8 + c)
		.attr("y2", 20)
    	.attr("transform", tf_element)
    	.style("stroke", color_intron)
    	.style("stroke-width", 2)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");	
	var position_t = off_set * 8 + c + 72;									
	s_l.append("path")
		.attr("d", s_sy)
		.attr("fill", "none")
		.attr("stroke","black")
		.attr("stroke-width", "1px")
		.style("opacity", "0.0")
		.attr("transform", function () {
							 return "translate(" + position_t + ",6)" + "rotate(90)";  })
		.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
	s_l.append("path")
		.attr("d", s_sy)
		.attr("fill", "none")
		.attr("stroke","black")
		.attr("stroke-width", "1px")
		.style("opacity", "0.0")
		.attr("transform", function () {
							 return "translate(" + position_t + ",22)" + "rotate(90)";  })
		.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
	s_l.append("text")
    	.attr("x", off_set * 9)
    	.attr("y", 10)
    	.style("font-size", "12px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.style("opacity", "0.0")
    	.attr("transform", tf_text)
    	.text("Type 5'").transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0"); 
	    
	//type 5
	s_l.append("line")
    	.attr("x1", off_set * 8 + c)
    	.attr("y1", off_set_y + 5)
		.attr("x2", off_set * 8 + c)
		.attr("y2", off_set_y + 21)
    	.attr("transform", tf_element)
    	.style("stroke", color_intron)
    	.style("stroke-width", 2)
    	.style("opacity", "0.0")
    	.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");	
	var position_f = off_set * 8 + c + 68;
	s_l.append("path")
		.attr("d", s_sy)
		.attr("fill", "none")
		.attr("stroke","black")
		.attr("stroke-width", "1px")
		.style("opacity", "0.0")
		.attr("transform", function () {
							 return "translate(" + position_f + "," + (off_set_y + 7) + ")" + "rotate(-90)";  })
		.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
	s_l.append("path")
		.attr("d", s_sy)
		.attr("fill", "none")
		.attr("stroke","black")
		.attr("stroke-width", "1px")
		.style("opacity", "0.0")
		.attr("transform", function () {
							 return "translate(" + position_f + "," + (off_set_y + 22) + ")" + "rotate(-90)";  })
		.transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
	s_l.append("text")
    	.attr("x", off_set * 9)
    	.attr("y", off_set_y + 10)
    	.style("font-size", "12px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.style("opacity", "0.0")
    	.attr("transform", tf_text)
    	.text("Type 3'").transition()
	    .delay(1500)
	    .duration(750)
	    .style("opacity", "1.0");
}

/* RM_ELEMENT_INFO_BOX
 * 
 * Rimuove gli elementi presenti nell'info_box al cambio
 * del gene. Gli elementi prima di essere eliminati,
 * svaniscono per effetto del decremento della trasparenza.
 */
function remove_element_info_box(){
	
	d3.select("#regions_selected")
		.transition()
    	.duration(750)
    	.style("opacity", "0.0")
		.remove(); 
                        
    d3.selectAll("#exon")
    	.transition()
        .duration(750)
        .style("opacity","1.0")
    	.attr("pointer-events", "yes")
      	.style("stroke-width", 0)
      	.style("fill", function() { return d3.rgb("#228B22"); });
                                   
    d3.selectAll("#exon_stripe")
    	.transition()
        .duration(750)
        .style("opacity","1.0")
    	.attr("pointer-events", "yes")
        .style("stroke-width", 0)
        .style("fill", 'url(#diagonalHatch)');
    d3.select("#exon_stripe_over")
    	.transition()
        .duration(750)
        .style("opacity", "0.0")
    	.remove();
                        
    d3.selectAll("#exon")
    	.transition()
        .duration(750)
        .style("opacity","1.0")
    	.attr("pointer-events", "yes")
        .style("stroke-width", 0)
        .style("fill", function() { return d3.rgb("#228B22"); });
                        
    d3.selectAll("#intron")
    	.transition()
        .duration(750)
        .style("opacity","1.0")
    	.style("stroke", "black")
    	.attr("pointer-events", "yes");
    
    d3.select("#title_sequence_box")
    	.transition()
    	.duration(750)
    	.style("opacity", "0.0")
    	.remove();    
    d3.select("#sequence_ex")
    	.transition()
    	.duration(750)
    	.style("opacity", "0.0")
    	.remove();
    d3.select("#sequence_in")
    	.transition()
    	.duration(750)
    	.style("opacity", "0.0")
    	.remove();
    d3.select("#table_title")
    	.transition()
    	.duration(750)
    	.style("opacity", "0.0")
    	.remove();
    d3.select("#table_start")
    	.transition()
    	.duration(750)
    	.style("opacity", "0.0")
    	.remove();
    d3.select("#table_end")
    	.transition()
    	.duration(750)
    	.style("opacity", "0.0")
    	.remove(); 
    	
    d3.select("#expande_info")
    	.style("border-color", function() { return "rgb(189, 195, 199)"; });
	d3.select("#isoform")
		.style("border-color", function() { return "rgba(34, 139, 34, 0.7)"; });    
}


/* SVG_INFO_BOX
 * 
 * Crea una finestra dove saranno visualizzati gli elementi
 * appartenenti alla selezione. Sono presenti i tasti per
 * fare un clean della finestra e riattivare la struttura del gene
 * e per visualizzare le informazione di start/end degli elementi 
 * visualizzati.
 */
function svg_info_box(){
	
	//dimensioni fisse della finestra degli elementi selezionati
    s_w = 750;
    s_h = 250;
    //la finestra viene centrata in base alla sua dimensione e alla
    //dimensione della finestra del browser
    var pos_w = ((window.innerWidth - s_w)/2);
    //vettore per il posizionamento            
    var p_s = {
    	pos : "absolute", 
    	left : pos_w + "px",
    	right : "10px",
    	top: "380px", 
    	bottom : "10px"
    };
                                   
    var s_i = set_svg("expande_info", s_w, s_h, p_s);
    s_i.attr("viewbox", function() { return "0 0" + s_w + s_h; });
    
    //bottone per cancellare il contenuto della finestra e riattivare
    //la struttura del gene    
    d3.select("body").append("button")
        .attr("id", "clear_button")
        .text("Clear")
        .style("top", "380px")
        .style("left", function () { return (pos_w + s_w + 7) + "px"; })
        .style("position", "absolute")
        .on("click", remove_element_info_box); 
    
    d3.select("body").append("button")
        .attr("id", "show_button")
        .text("Show Info")
        .style("top", "410px")
        .style("left", function () { return (pos_w + s_w + 7) + "px"; })
        .style("position", "absolute")
        .on("click", function() {
            d3.select("#table_title").attr("visibility", "visible")
            	.style("opacity", "0.0")
            	.transition()
            	.duration(750)
            	.style("opacity", "1.0");
            d3.select("#table_start").attr("visibility", "visible")
            	.style("opacity", "0.0")
            	.transition()
            	.duration(750)
            	.style("opacity", "1.0");
            d3.select("#table_end").attr("visibility", "visible")
            	.style("opacity", "0.0")
            	.transition()
            	.duration(750)
            	.style("opacity", "1.0");
        }); 
    
    return s_i;
}


/* ELEMENT_SELECTED_EXON
 * e_i -> elemento "esone" selezionato
 * s -> posizione di "start" dell'elemento
 * e -> posizione di "end" dell'elemento
 * 
 * Disegna un nuovo esone nella stessa posizione di quello 
 * selezionato. Si ottiene così l'effetto di evidenziare 
 * solo l'esone selezionato
 */
function element_selected_exon(e_i, s, e){
    
    //traslazione elemento
    var tf_e = d3.svg.transform()
		.translate(function (d) { return [s, 0]; });
    
    var ex = e_i.append("rect")
        .attr("id", "element_info_exon")
        .attr("x", s)
        .attr("y", 0)
        .attr("width", function() { return e - s; })
        .attr("height", 75)
        .style("fill", function() { return d3.rgb("#228B22").brighter(3); });
    return ex;
}


/* CONNECT_EXON
* exon -> esone sulla struttura del gene
* sequence -> sequenza del'esone
*
* Crea un link tra l'esone selezionato e la corrispondente sequenza nucleotidica
*/
function connect_exon(exon, sequence){
		
	var xr_s = +exon.attr("x");
	var yr_s = +exon.attr("y") + +exon.attr("height") + 1;
	var xt_s = +sequence.attr("x") + 4;
	var yt_s = +sequence.attr("y") - 7;
	
	var xr_e = +exon.attr("x") + +exon.attr("width");
	var yr_e = +exon.attr("y") + +exon.attr("height") + 1;
	var xt_e = +sequence.attr("x") + sequence.text().length*10 + 2;
	var yt_e = +sequence.attr("y") - 7;

    point = [{ "source" : { "x" : xr_s + margin_isoform.left, "y" : yr_s}, 
    		   "target" : { "x" : xt_s + (margin_isoform.left * 20), "y" : yt_s + 135}},
    		 { "source" : { "x" : xr_e + margin_isoform.left, "y" : yr_e}, 
    		   "target" : { "x" : xt_e + (margin_isoform.left * 20), "y" : yt_e + 135}} ];
    		  
    var diagonal = d3.svg.diagonal()
    	.source(function(d) { return {"x":d.source.y, "y":d.source.x}; })            
    	.target(function(d) { return {"x":d.target.y, "y":d.target.x}; })
        .projection(function(d) { return [+d.y, +d.x + 50]; });
        	
    svg_box.selectAll(".link")
    	.data(point)
    	.enter().append("path")
    	.attr("class", "link")
        .attr("d", diagonal);      
}


/* CONNECT_INTRON
* intron -> introne sulla struttura del gene
* sequence -> sequenza del'esone
*
* Crea un link tra l'introne selezionato e la corrispondente sequenza nucleotidica
*/
function connect_intron(intron, sequence){
	
	var xr_s = +intron.attr("x1");
	var yr_s = +intron.attr("y1");
	var xt_s = +sequence.attr("x") + 4;
	var yt_s = +sequence.attr("y") - 7;
	
	var xr_e = +intron.attr("x2");
	var yr_e = +intron.attr("y2");
	var xt_e = +sequence.attr("x") + sequence.text().length*10;
	var yt_e = +sequence.attr("y") - 7;
	
	var off_set_s = off_set_ex - sequence.text().length*12 + 2;
	var off_set_e = off_set_ex - sequence.text().length*10;

    point = [{ "source" : { "x" : xr_s + margin_isoform.left, "y" : yr_s}, 
    		   "target" : { "x" : xt_s + (margin_isoform.left * 20) + off_set_s, "y" : yt_s + 135}},
    		 { "source" : { "x" : xr_e + margin_isoform.left, "y" : yr_e}, 
    		   "target" : { "x" : xt_e + (margin_isoform.left * 20) + off_set_e, "y" : yt_e + 135}} ];
    		  
    var diagonal = d3.svg.diagonal()
    	.source(function(d) { return {"x":d.source.y, "y":d.source.x}; })            
    	.target(function(d) { return {"x":d.target.y, "y":d.target.x}; })
        .projection(function(d) { return [+d.y, +d.x + 50]; });
        	
    svg_box.selectAll(".link")
    	.data(point)
    	.enter().append("path")
    	.attr("class", "link")
        .attr("d", diagonal);      
}


/* ELEMENT_SELECTED_INTRON
 * s -> posizione di "start" dell'elemento
 * e -> posizione di "end" dell'elemento
 * 
 * Disegna un nuovo introne nella stessa posizione di quello 
 * selezionato. Si ottiene così l'effetto di evidenziare 
 * solo l'introne selezionato
 */
function element_selected_intron(s, e){
     
    var intron_selected = d3.select("#introns").append("line")
        .attr("id", "element_info_intron")
		.attr("x1", s)
		.attr("y1", 35)
		.attr("x2", e)
		.attr("y2", 35)
		.style("stroke", "black")
		.style("stroke-width", 6);
		
	return intron_selected;
}


function display_info(s_i, domain, elements, r, x){
    
    //elementi estratti dalla selezione
    var exons_info = elements[0];
    var introns_info = elements[1];
    console.log(exons_info);
    console.log(introns_info);
    
    var color_exon = function() { return d3.rgb("#228B22"); };
    var color_intron = function() { return d3.rgb("black"); };
    
        
    tipElement = d3.tip()
        .attr('class', 'd3-tip')
        .offset([10, 0])
        .direction("e")
        .html(function(d) {
                return "<br><strong>Start:</strong> <span style='color:yellow'>" + 
                        d.start + "</span>" +
                        "<br><strong>End:</strong> <span style='color:yellow'>" + 
                        d.end + "</span>";
                    
        });
      
    //variabili per le operazioni di trasformazione 
    var transf = {
    	
    	//traslazione esoni   
   		t_e : d3.svg.transform()
        	.translate(function (d, i) { return [30, i * 45]; })
        	.scale(function (d, i) { return [2.8, 1]; }),
        //traslazione introni
    	t_i : d3.svg.transform()
        	.translate(function (d, i) { return [30, (i * 20) + (exons_info.length * 45)]; })
        	.scale(function () { return [2, 1]; }),
        //traslazione contenitore elementi 
    	tf_g : d3.svg.transform()
        	.translate(function (d, i) { return [15, 20]; }),
        tf_table_title : d3.svg.transform()
        	.translate(function (d, i) { return [s_w*(1/2), 0]; }),
        tf_table_start : d3.svg.transform()
        	.translate(function (d, i) { return [s_w*(1/2), 45]; }),
        tf_table_end : d3.svg.transform()
        	.translate(function (d, i) { return [s_w*(3/4), 45]; })
       };
       
    var table_text = [];
    for(k = 0; k < exons_info.length; k++)
    	table_text.push(exons_info[k]);
    if(introns_info != null)
    	for(k = 0; k < introns_info.length; k++)
    		table_text.push(introns_info[k]);
       
    var table_title = s_i.append("g")
        .attr("id", "table_title")
        .attr("transform", transf.tf_table_title)
        .attr("visibility", "hidden");
    
    table_title.append("text")
    	.attr("x", s_w*(1/8))
    	.attr("y", 15)
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "blue")
    	.text("Start");
    table_title.append("text")
    	.attr("x", s_w*(1/3) + 30)
    	.attr("y", 15)
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "blue")
    	.text("End");
    table_title.append("line")
    	.attr("x1", (s_w*(1/8)) - 30)
    	.attr("y1", 20)
    	.attr("x2", (s_w*(1/8)) + 250)
    	.attr("y2", 20)
    	.style("stroke", "black")
    	.style("stroke-width", "1px");
    	
    var table_start = s_i.append("g")
        .attr("id", "table_start")
        .attr("transform", transf.tf_table_start)
        .attr("visibility", "hidden");
    
    table_start.selectAll("text")
    	.data(table_text)
    	.enter().append("text")
    	.attr("id", function(d) { return "text_" + d.id; })
    	.attr("transform", function(d, i) { 
    						console.log(d.start.toString().length);
    						if(d.pattern == null)
    							return "translate(" + (s_w*(1/8) - (d.start.toString().length)) + "," +  i * 45 + ")";
    						else
    							return "translate(" + (s_w*(1/8) - (d.start.toString().length)) + "," + ((i * 20) + (exons_info.length * 35)) + ")"; })
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.text(function(d) { return d.start; });
    
    var table_end = s_i.append("g")
        .attr("id", "table_end")
        .attr("transform", transf.tf_table_end)
        .attr("visibility", "hidden");
    
    table_end.selectAll("text")
    	.data(table_text)
    	.enter().append("text")
    	.attr("id", function(d) { return "text_" + d.id; })
    	.attr("transform", function(d, i) { 
    						if(d.pattern == null)
    							return "translate(" + (s_w*(1/8) - (d.start.toString().length)) + "," +  i * 45 + ")";
    						else
    							return "translate(" + (s_w*(1/8) - (d.start.toString().length)) + "," + ((i * 20) + (exons_info.length * 35)) + ")"; })
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.text(function(d) { return d.end; });
    
    var g = s_i.append("g")
        .attr("id", "regions_selected")
        .attr("transform", transf.tf_g);        
    g.selectAll("rect")
        .data(exons_info)
        .enter().append("rect")
        .attr("width", function(d) { return domain(d.end) - domain(d.start); })
        .attr("height", 40)
        .style("fill", color_exon)
        .style("opacity", "0.0")
        .attr("transform", transf.t_e)
        .on("mouseover", function(d) { 
                            d3.select(this).style('cursor', 'crosshair');
                            //esone selezionato nella struttura
                            var s = element_selected_exon(r, x(d.start), x(d.end));
                            
                            seq_id = "#sequence_ex_" + d.id;
                            //sequenza nucleotidica corrispondente
                            var t = d3.select(seq_id);
                            t.style("fill", "#822222");
                            connect_exon(s, t);
                            d3.selectAll("#text_" + d.id)
                            	.style("fill", "#822222");
                            })
        .on("mouseout", function(d) { 
                            d3.select(this).style('cursor', 'default');
                            d3.select("#element_info_exon").remove();
                            
                            seq_id = "#sequence_ex_" + d.id;
                            d3.select(seq_id).style("fill", "black");
                            d3.selectAll(".link").remove();
                            d3.selectAll("#text_" + d.id)
                            	.style("fill", "black");
                           })
        .transition()
        .duration(750)
        .style("opacity","1.0");
           
    if(introns_info != null){
    	g.selectAll("line")
    		.data(introns_info)
			.enter().append("line")
			.attr("x1", function(d) { return domain(d.start); })
			.attr("y1", 35)
			.attr("x2", function(d) { return domain(d.end); })
			.attr("y2", 35)
			.attr("transform", transf.t_i)
			.style("stroke", color_intron)
			.style("stroke-width", 8)
			.style("opacity", "0.0")
            .on("mouseover", function(d, i) { 
                                d3.select(this).style('cursor', 'crosshair');
                                //traslazione testo pattern
                                var tf_info_text = d3.svg.transform()
                                    .translate(function () { 
                                        return [30, (i * 20) + (exons_info.length * 45)]; })
                                    .scale(function () { return [2, 1]; });
                                    
                                g.append("text")
                                    .attr("x", domain(d.start) - 15)
                                    .attr("y", 38)
                                    .style("font-size", "10px")
                                    .attr("transform", tf_info_text)
                                    .style("font-family", "Arial, Helvetica, sans-serif")
                                    .style("fill", "black")
                                    .text(d.pattern.slice(0,2).toUpperCase());
                                g.append("text")
                                    .attr("x", domain(d.end + 10))
                                    .attr("y", 38)
                                    .style("font-size", "10px")
                                    .style("font-family", "Arial, Helvetica, sans-serif")
                                    .attr("transform", tf_info_text)
                                    .style("fill", "black")
                                    .text(d.pattern.slice(2,4).toUpperCase());
                                    
                                var s = element_selected_intron(x(d.start), x(d.end));
                                seq_id = "#sequence_in_" + d.id;
                                var t = d3.select(seq_id);
                                t.style("fill", "#822222");
                                connect_intron(s, t);
                                d3.selectAll("#text_" + d.id)
                            	.style("fill", "#822222"); })
            .on("mouseout", function(d) { 
                                d3.select(this).style('cursor', 'default');
                                d3.select("#element_info_intron").remove();
                                seq_id = "#sequence_in_" + d.id;
                                d3.select(seq_id).style("fill", "black");
                                g.selectAll("text").remove();
                                d3.selectAll(".link").remove();
                                d3.selectAll("#text_" + d.id)
                            	.style("fill", "black"); })
            .transition()
        	.duration(750)
        	.style("opacity","1.0");
	}   	
}


/* DISPLAY_INFO_STRIPE
 * s_i -> finestra creata per visualizzare le informazioni
 * domain -> dominio della finestra di visualizzazione degli elementi
 * elements -> elementi estratti dalla struttura del gene
 *             in base alla selezione.
 * r -> contenitore degli esoni
 * x -> variabile che contiene la funzione per il range della finestra
 * 
 * Visualizza gli elementi appartenenti alla selezione 
 * (per esone 'conservative').
 */     
function display_info_stripe(s_i, domain, elements, r, x){
    
    //elementi estratti dalla selezione
    var exons_info = elements[0];
    var introns_info = elements[1];
    
    
    var color_exon = function() { return d3.rgb("#228B22"); };
    var color_intron = function() { return d3.rgb("black"); };
    
    var transf = {
    	
    	//traslazione esoni   
   		t_e : d3.svg.transform()
        	.translate(function (d, i) { return [30, i * 45]; })
        	.scale(function (d, i) { return [2.8, 1]; }),
        //traslazione introni
    	t_i : d3.svg.transform()
        	.translate(function (d, i) { return [30, (i * 20) + (exons_info.length * 45)]; })
        	.scale(function () { return [2, 1]; }),
        //traslazione contenitore elementi 
    	tf_g : d3.svg.transform()
        	.translate(function (d, i) { return [15, 20]; }),
        tf_table_title : d3.svg.transform()
        	.translate(function (d, i) { return [s_w*(1/2), 0]; }),
        tf_table_start : d3.svg.transform()
        	.translate(function (d, i) { return [s_w*(1/2), 45]; }),
        tf_table_end : d3.svg.transform()
        	.translate(function (d, i) { return [s_w*(3/4), 45]; })
       };
    
    var table_text = [];
    for(k = 0; k < exons_info.length; k++)
    	table_text.push(exons_info[k]);
    if(introns_info != null)
    	for(k = 0; k < introns_info.length; k++)
    		table_text.push(introns_info[k]);
    
    //head tabella   
    var table_title = s_i.append("g")
        .attr("id", "table_title")
        .attr("transform", transf.tf_table_title)
        .attr("visibility", "hidden");
    table_title.append("text")
    	.attr("x", s_w*(1/8))
    	.attr("y", 15)
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "blue")
    	.text("Start");
    table_title.append("text")
    	.attr("x", s_w*(1/3) + 30)
    	.attr("y", 15)
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "blue")
    	.text("End");
    table_title.append("line")
    	.attr("x1", (s_w*(1/8)) - 30)
    	.attr("y1", 20)
    	.attr("x2", (s_w*(1/8)) + 250)
    	.attr("y2", 20)
    	.style("stroke", "black")
    	.style("stroke-width", "1px");
    
    //contenuto tabella	
    var table_start = s_i.append("g")
        .attr("id", "table_start")
        .attr("transform", transf.tf_table_start)
        .attr("visibility", "hidden"); 
    table_start.selectAll("text")
    	.data(table_text)
    	.enter().append("text")
    	.attr("id", function(d) { return "text_" + d.id; })
    	.attr("transform", function(d, i) { 
    						if(d.pattern == null)
    							return "translate(" + (s_w*(1/8) - 10) + "," +  i * 45 + ")";
    						else
    							return "translate(" + (s_w*(1/8) - 10) + "," + ((i * 20) + (exons_info.length * 35)) + ")"; })
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.text(function(d) { return d.start; });  
    var table_end = s_i.append("g")
        .attr("id", "table_end")
        .attr("transform", transf.tf_table_end)
        .attr("visibility", "hidden");   
    table_end.selectAll("text")
    	.data(table_text)
    	.enter().append("text")
    	.attr("id", function(d) { return "text_" + d.id; })
    	.attr("transform", function(d, i) { 
    						if(d.pattern == null)
    							return "translate(" + (s_w*(1/8) - 10) + "," +  i * 45 + ")";
    						else
    							return "translate(" + (s_w*(1/8) - 10) + "," + ((i * 20) + (exons_info.length * 35)) + ")"; })
    	.style("font-size", "16px")
    	.style("font-family", "Arial, Helvetica, sans-serif")
    	.style("fill", "black")
    	.text(function(d) { return d.end; });
    
    //esoni
    var g = s_i.append("g")
        .attr("id", "regions_selected")
        .attr("transform", transf.tf_g);        
    g.selectAll("rect")
        .data(exons_info)
        .enter().append("rect")
        .attr("width", function(d) { return domain(d.end) - domain(d.start); })
        .attr("height", 40)
        .style("fill", color_exon)
        .style("opacity", "0.0")
        .attr("transform", transf.t_e)
        .on("mouseover", function(d) { 
        					d3.select(this).style('cursor', 'crosshair');
                            var s = element_selected_exon(r, x(d.start), x(d.end)); 
                            
                            seq_id = "#sequence_ex_" + d.id;
                            var t = d3.select(seq_id);
                            t.style("fill", "#822222");
                            connect_exon(s, t);
                            d3.selectAll("#text_" + d.id)
                            	.style("fill", "#822222");})
        .on("mouseout", function(d) { 
        					d3.select(this).style('cursor', 'default');
        					d3.select("#element_info_exon").remove();
        					
        					seq_id = "#sequence_ex_" + d.id;
                            d3.select(seq_id).style("fill", "black");
                            d3.selectAll(".link").remove();
                            d3.selectAll("#text_" + d.id)
                            	.style("fill", "black"); })
        .transition()
        .duration(750)
        .style("opacity","1.0");
    
    //introni    
    if(introns_info != null){
    	g.selectAll("line")
    		.data(introns_info)
			.enter().append("line")
			.attr("x1", function(d) { return domain(d.start); })
			.attr("y1", 35)
			.attr("x2", function(d) { return domain(d.end); })
			.attr("y2", 35)
			.attr("transform", transf.t_i)
			.style("stroke", color_intron)
			.style("stroke-width", 8)
			.style("opacity", "0.0")
			.on("mouseover", function(d, i) { 
                                d3.select(this).style('cursor', 'crosshair');
                                var s = element_selected_intron(x(d.start), x(d.end));
                                //traslazione testo pattern
                                var tf_info_text = d3.svg.transform()
                                    .translate(function () { 
                                        return [30, (i * 20) + (exons_info.length * 45)]; })
                                    .scale(function () { return [2, 1]; });
                                    
                                g.append("text")
                                    .attr("x", domain(d.start) - 15)
                                    .attr("y", 38)
                                    .attr("font-size", "10px")
                                    .attr("transform", tf_info_text)
                                    .style("fill", "black")
                                    .text(d.pattern.slice(0,2).toUpperCase());
                                g.append("text")
                                    .attr("x", domain(d.end + 10))
                                    .attr("y", 38)
                                    .attr("font-size", "10px")
                                    .attr("transform", tf_info_text)
                                    .style("fill", "black")
                                    .text(d.pattern.slice(2,4).toUpperCase());
                                    
                                seq_id = "#sequence_in_" + d.id;
                                var t = d3.select(seq_id);
                                t.style("fill", "#822222");
                                connect_intron(s, t);
                                d3.selectAll("#text_" + d.id)
                            	.style("fill", "#822222"); })
            .on("mouseout", function(d) { 
                                d3.select(this).style('cursor', 'default');
                                d3.select("#element_info_intron").remove();
                                seq_id = "#sequence_in_" + d.id;
                                d3.select(seq_id).style("fill", "black");
                                g.selectAll("text").remove();
                                d3.selectAll(".link").remove();
                                d3.selectAll("#text_" + d.id)
                            	.style("fill", "black"); })
            .transition()
        	.duration(750)
        	.style("opacity","1.0");   
	}     
}


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
 * alle coordinate della selezione. Restituisce una struttura di array.
 */
function check_structure_element(regions_info, c){
    
    //elementi che verranno estratti dalla selezione 
    var element = {
    	r_i : [],
    	i_i : [],
    	c_i : []
    };
    
    for(rg = 0; rg < regions_info.length; rg++)
        element.r_i.push(regions[regions_info[rg]]);
    
    for(i = 0; i < introns_restruct.length; i++)
        if((c >= introns_restruct[i].start) & (c <= introns_restruct[i].end))
            element.i_i.push(introns_restruct[i]);
            
    if(element.r_i.length != 0)
        element.c_i.push(element.r_i);
    if(element.i_i.length != 0)
        element.c_i.push(element.i_i);
    
    return element.c_i;   
}


/* DRAW_EXONS
 * box -> variabile che contiente l'elemento "svg"
 * exons -> struttura dati degli esoni
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Disegna gli esoni differenziandoli in base al campo "alternative". In base
 * alla posizione del mouse e all'evento connesso ad esso, richiama le funzioni 
 * per la visualizzazione dell'elemento selezionato.
 */
function draw_exons(box, exons, x_scale){
	
	//array per gli esoni "conservative"
	var exons_stripe = [];
	var exons_h = 75;
	
	//colori per gli elementi
	var color_exon = function() { return d3.rgb("#228B22"); };
	var color_exon_after = function() { return d3.rgb("#808080"); };
    var color_intron = function() { return d3.rgb("black"); };
    var color_intron_after = function() { return d3.rgb("#808080"); };
	
	for(k = 0; k < exons.length; k++)
	   if(exons[k].alternative == false)
	       exons_stripe.push(exons[k]);
    	  			  	
	//variabile per traslare gli esoni
	var tf = d3.svg.transform()
		.translate(function (d) { return [x_scale(d.start), 0]; });
	
	//contenitore degli esoni
	var rect_exons = box.append("g")
		.attr("id", "exons")
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top/2 + ")");
	//aggiunge i blocchi "esoni"
	rect_exons.selectAll("rect")
		.data(exons)
		.enter()
		.append("rect")
		.attr("id", "exon")
		.attr("width", function(d) { return x_scale(d.end) - x_scale(d.start) - 1; })
		.attr("height", "0")
		.attr("rx", 3)
		.attr("ry", 3)
		.style("fill", color_exon)
		.style("opacity", 1.0)
		.attr("transform",tf)
		.on("click", function(d){ 
		               if(d.alternative == true){
		                  d3.selectAll("#exon")
						    .style("fill", color_exon_after)
							.attr("pointer-events", "none");
						  
						  d3.selectAll("#exon_stripe")
							.attr("pointer-events", "none");
							
						  d3.selectAll("#intron")
                            .style("stroke", color_intron_after)
                            .attr("pointer-events", "none");
                          
                          d3.select(this)
                          .style("fill", color_exon);
									
						  //coordinate della posizione del mouse al momento del "click"
						  var coord_x = x_scale.invert(d3.event.pageX);					  
						  if((coord_x > x_scale(d.start)) | (coord_x < x_scale(d.end))){
						      info_structure = check_structure_element(d.regions, coord_x);  	
						      x_info = window_info_scale(info_structure, s_h);
						      display_info(svg_info, x_info, info_structure, rect_exons, x_scale);
						      sequence_box(box, info_structure);
						      d3.select("#expande_info")
						      	.style("border-color", function() { return "rgba(34, 139, 34, 0.7)"; });
						      d3.select("#isoform")
						      	.style("border-color", function() { return "rgb(189, 195, 199)"; });
						      }
        			      }})							  		
		.on("mouseover", function() { d3.select(this).style('cursor', 'cell'); })
		.on("mouseout", function() { d3.select(this).style('cursor', 'default'); })
	    .transition()
	    .duration(750)
	    .attr("height", exons_h);
	
	//esoni alternative = false								 
	var rect_exons_stripe = box.append("g")
        .attr("id", "exons")
        .attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top/2 + ")");
    //aggiunge i blocchi "esoni alternative"
    rect_exons_stripe.selectAll("rect")
        .data(exons_stripe)
        .enter()
        .append("rect")
        .attr("id", "exon_stripe")
        .attr("width", function(d) { return x_scale(d.end) - x_scale(d.start) - 1; })
        .attr("height", "0")
        .attr("rx", 3)
		.attr("ry", 3)
        .style("fill", 'url(#diagonalHatch)')
        .attr("transform",tf) 
        .on("click", function(d){ 
                        d3.selectAll("#exon_stripe")
                          .style("fill", color_exon_after)
                          .attr("pointer-events", "none");
                        d3.selectAll("#exon")
                          .style("fill", color_exon_after)
                          .attr("pointer-events", "none");
                        d3.selectAll("#intron")
                            .style("stroke", color_intron_after)
                            .attr("pointer-events", "none");
                                    
                        d3.select(this)
                          .style("fill", color_exon);
                          
                        rect_exons_stripe.append("rect")
                        	.attr("id", "exon_stripe_over")
        					.attr("width", x_scale(d.end) - x_scale(d.start) - 1)
        					.attr("height", 75)
        					.attr("rx", 3)
							.attr("ry", 3)
        					.style("fill", 'url(#diagonalHatch)')
        					.attr("transform", "translate(" + x_scale(d.start) + ",0)");
        		        
        		        //coordinate della posizione del mouse al momento del "click"       
                        var coord_x = x_scale.invert(d3.event.pageX);
                        if((coord_x > x_scale(d.start)) | (coord_x < x_scale(d.end))){
                        	info_structure = check_structure_element(d.regions, coord_x);  
                            x_info = window_info_scale(info_structure, s_h);
                            display_info_stripe(svg_info, x_info, info_structure, rect_exons_stripe, x_scale);
                            sequence_box(box, info_structure);
                            d3.select("#expande_info")
						      	.style("border-color", function() { return "rgba(34, 139, 34, 0.7)"; });
						    d3.select("#isoform")
						      	.style("border-color", function() { return "rgb(189, 195, 199)"; });
						    }
                    })                           
        .on("mouseover", function() { d3.select(this).style('cursor', 'cell'); })
        .on("mouseout", function() { d3.select(this).style('cursor', 'default'); })
        .transition()
        .duration(750)
        .attr("height", exons_h);
        		
	return rect_exons;	
}



/* DRAW_INTRONS
 * box -> variabile che contiente l'elemento "svg"
 * introns -> struttura dati degli introni
 * x_scale -> variabile che contiente la funzione per il range e il dominio di
 * 			  visualizzazione
 * 
 * Disegna gli introni.
 */
function draw_introns(box, introns, x_scale){
    
    var color_intron = function() { return d3.rgb("black"); };
    
    //contenitore degli introni
	var line_introns = box.append("g")
		.attr("id", "introns")
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top/2 + ")");
		
	line_introns.selectAll("line")
		.data(introns)
		.enter().append("line")
		.attr("id", "intron")
		.attr("x1", function(d) { return x(d.start); })
		.attr("y1", 35)
		.attr("x2", function(d) { return x(d.start); })
		.attr("y2", 35)
		.style("stroke", color_intron)
		.style("stroke-width", 6)
		.transition()
		.delay(750)
		.duration(750)
		.attr("x2", function(d) { return x(d.end); });
			
	return line_introns;
}



/* CLONE_TRIANGLE_UP (DEPRECATED)
 * svg -> variabile che contiene l'elemento "svg"
 * obj -> oggetto da clonare
 * 
 * Clona l'oggetto contenuto in "obj" e lo aggiunge alla finestra di 
 * visualizzazione contenuta nella variabile "svg".
 * Il contenitore "g" del "segnale alto" della tipologia degli splice sites 
 * viene clonato e riutilizzato per aggiungere un "segnale basso".
 */
function clone_svg_element(svg, obj) {
	
	var triangle_down = svg.append("use")
    	.attr("xlink:xlink:href","#" + obj.attr("id"));
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
	
	var color_s_s = function() { return d3.rgb("black"); };
	
	//variabile per i simboli della tipologia 
	//di splice_sites
	var s_sy = d3.svg.symbol()
		.type('triangle-up')
		.size(20);
					
	//contenitore degli splice sites
	var splice_sites = box.append("g")
		.attr("id", "splice_sites")
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top/2 + ")");
		
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
		.style("stroke", color_s_s)
		.style("strole-width", "2px")
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
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top/2 + ")");
										
	triangle_up.selectAll("path")
		.data(s_s)
		.enter().append("path")
		.attr("id", function(d, i) { return "up" + i; })
		.attr("d", s_sy)
		.attr("fill", "none")
		.attr("stroke","black")
		.attr("stroke-width", "1px")
		.style("opacity", "0.0")
		.attr("transform", function (d) {
		                      if((d.type == 3) | (d.type == "term"))
							     return "translate(" + (x_scale(d.position) - 3) + ",-27)" + "rotate(-90)";
							  else
							     if((d.type != "unknow") | (d.type == "init"))
								    return "translate(" + (x_scale(d.position) + 3) + ",-27)" + "rotate(90)"; })
	    .transition()
        .delay(1500)
        .duration(750)
        .style("opacity", "1.0");
    
    var triangle_down = box.append("g")
		.attr("id", "triangle_down")
		.attr("transform", "translate(" + margin_isoform.left + "," + margin_isoform.top/2 + ")");
											
	triangle_down.selectAll("path")
		.data(s_s)
		.enter().append("path")
		.attr("id", function(d, i) { return "up" + i; })
		.attr("d", s_sy)
		.attr("fill", "none")
		.attr("stroke","black")
		.attr("stroke-width", "1px")
		.style("opacity", "0.0")
		.attr("transform", function (d) {
		                      if((d.type == 3) | (d.type == "term"))
							     return "translate(" + (x_scale(d.position) - 3) + ", 106)" + "rotate(-90)";
							  else
							     if((d.type != "unknow") | (d.type == "init"))
								    return "translate(" + (x_scale(d.position) + 3) + ", 106)" + "rotate(90)"; })
	    .transition()
        .delay(1500)
        .duration(750)
        .style("opacity", "1.0");
        
	/*//viene clonato triangle_up (DEPRECATED)									 
	var td = clone_svg_element(box, triangle_up);
	td.attr("transform", "translate(0, 136)")
    	.attr("id", "triangle_down")
    	.attr("opacity", "0.0")
    	.transition()
        .delay(1500)
        .duration(750)
        .style("opacity", "1.0");*/
			
	return splice_sites;
}


/* SEQUENCE_BOX
 * isoform_box -> contenitore della struttura del gene
 * seq_info -> array contenente gli elementi
 * 
 * Visualizza le sequenze nucleotidiche degli elementi selezionati
 */
function sequence_box(isoform_box, seq_info){
	
	var color_title = function() { return d3.rgb("blue"); };
	var color_sequence = function() { return d3.rgb("black"); };
	
	//array delle posizioni in base al numero di sequenze 
	//da visualizzare
	var position = {
		y_pos : 190,
		x_pos_ex : 250,
		x_pos_exx : 350,
		x_pos_in : 150
	};
    
    isoform_box.append("text")
        .attr("id", "title_sequence_box")
        .attr("x", 0)
        .attr("y", 30)
        .style("font-family", "Arial, Helvetica, sans-serif")
        .style("font-size", "20px")
        .style("fill", color_title)
        .text("Nucleic sequence:")
        .attr("transform", "translate(" + margin_isoform.left + "," + position.y_pos + ")")
        .style("opacity", "0.0")
        .transition()
        .duration(750)
        .style("opacity","1.0");
    
    //offset tra sequenze di esoni e introni 
    if(seq_info[0].length > 1)   	
		off_set_ex = seq_info[0].length * position.x_pos_ex;
	else
		off_set_ex = seq_info[0].length * position.x_pos_exx;
	if(seq_info[1] != null)
		var off_set_in = seq_info[1].length * position.x_pos_in;
	
	//box per le sequenze degli esoni
	var sequence_ex = isoform_box.append("g")
		.attr("id", "sequence_ex")
		.attr("transform", "translate(" + (margin_isoform.left * 20) + "," + position.y_pos + ")");
    
    //box per le sequenze degli introni
    var sequence_in = isoform_box.append("g")
        .attr("id", "sequence_in")
        .attr("transform", "translate(" + off_set_ex + "," + position.y_pos + ")");
   
    sequence_ex.selectAll("text")
    	.data(seq_info[0])
    	.enter().append("text")
    	.attr("id", function(d) { return "sequence_ex_" + d.id; })
        .attr("x", function (d, i) { 
        			var canvas = document.createElement('canvas');
					var ctx = canvas.getContext("2d");
					ctx.font = "20px Arial";        
					var width = ctx.measureText(d.sequence).width;
        			return (i * width); })
        .attr("y", 30)
        .style("font-size", "16px")
        .style("font-family", "Arial, Helvetica, sans-serif")
        .style("fill", color_sequence)
        .style("opacity", "0.0")
        .text(function(d) { return d.sequence.toUpperCase(); })
        .transition()
        .duration(750)
        .style("opacity","1.0");
    
    //solo se gli introni appartengono alla selezione    
    if(seq_info[1] != null)
        sequence_in.selectAll("text")
            .data(seq_info[1])
            .enter().append("text")
            .attr("id", function(d) { return "sequence_in_" + d.id; })
            .attr("x", function (d, i) {
            			var canvas = document.createElement('canvas');
						var ctx = canvas.getContext("2d");
						ctx.font = "23.5px Arial";        
						var width = ctx.measureText(d.prefix + d.suffix).width;
        				return (i * width); })
            .attr("y", 30)
            .style("font-size", "16px")
            .style("font-family", "Arial, Helvetica, sans-serif")
            .style("fill", color_sequence)
            .style("opacity", "0.0")
            .text(function(d) { return (d.prefix + d.suffix).toUpperCase(); })
            .transition()
        	.duration(750)
        	.style("opacity","1.0");        
}


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
    
    //rimozione di tutti i box
    var g = d3.select("#isoform").selectAll("g");
    g.remove();
    
    //rimozione del titolo
    d3.select("#title").remove();
    
    //pulisce la finestra degli elementi selezionati
    remove_element_info_box();
    
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
       .style("left", function () { return margin_isoform.left + "px"; })
       .style("position", "absolute");
    s_g.html(function() { return '<select><option value="ATP6AP1example1">ATP6AP1_ex1</option>' + 
    							   '<option value="ATP6AP1example2">ATP6AP1_ex2</option>' +
                                   '<option value="ATP6AP1example3">ATP6AP1_ex3</option>' + 
                                   '</select>'; }); 
    
    s = d3.select("select")
            .on("change", change_gene);
    
    s.property("value", "ATP6AP1example2");
    
    //texture esoni alternative
    pattern_exons();
    
    //dimensioni e posizione finestra della struttura
    //del gene
    var height_isoform = 250;
    var pos_box = {
    	pos: "absolute",
    	left : "10px",
    	right : "10px",
    	top : "110px",
    	bottom : "10px"
    };
    svg_box = set_svg("isoform", width - margin_isoform.right - margin_isoform.left, height_isoform, pos_box);
        
    init();   
    
    svg_info = svg_info_box();  
    
    legend_box();                         
}


/* DISPLAY_GENE
 * 
 * Visualizza il nome del gene di cui sarà visualizzata
 * la struttura.
 */
function display_gene(){
    
    var color_gene = function() { return d3.rgb("black"); };
    
    //variabile per il titolo
    var title = "";
    
    title = title.concat(original_info.gene + " gene structure");
    
    //array per il posizionamento
    var pos_title = {
    	pos : "absolute",
    	left : margin_isoform.left + "px",
    	right : margin_isoform.right + "px",
    	top : "40px",
    	bottom : "10px"
    };
    //altezza finestra
    var svg_height = 50;
    
    var svg_title = set_svg("title", width / 2, svg_height, pos_title);
    svg_title.append("text")
       .attr("id", "title")
       .attr("x", 10)
       .attr("y", 35)
       .style("font-family", "Arial, Helvetica, sans-serif")
       .style("font-size", "30px")
       .style("fill", color_gene)
       .style("opacity", "0.0")
       .text(title)
       .transition()
       .duration(1000)
       .style("opacity", "1.0");
      
	svg_title.append("text")
       .attr("x", 425)
       .attr("y", 15)
       .style("font-family", "Arial, Helvetica, sans-serif")
       .style("font-size", "12px")
       .style("fill", color_gene)
       .style("opacity", "0.0")
       .text(function () { 
       				var s = "sequence_id -> ";
       				return s + original_info.sequence_id; })
       .transition()
       .duration(1000)
       .style("opacity", "1.0");
       
    svg_title.append("text")
       .attr("x", 425)
       .attr("y", 30)
       .style("font-family", "Arial, Helvetica, sans-serif")
       .style("font-size", "12px")
       .style("fill", color_gene)
       .style("opacity", "0.0")
       .text(function () { 
       				var s = "program_version -> ";
       				return s + original_info.program_version; })
       .transition()
       .duration(1000)
       .style("opacity", "1.0");
    
    svg_title.append("text")
       .attr("x", 425)
       .attr("y", 45)
       .style("font-family", "Arial, Helvetica, sans-serif")
       .style("font-size", "12px")
       .style("fill", color_gene)
       .style("opacity", "0.0")
       .text(function () { 
       				var s = "file_format -> ";
       				return s + original_info.file_format_version; })
       .transition()
       .duration(1000)
       .style("opacity", "1.0");
}


/* COPY_REGIONS
 * s -> struttura da copiare
 * 
 * Copia il vettore delle regioni che poi sara modificato
 * per scalare le dimensioni degli elementi
 */
function copy_regions(s){
    
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

/* COPY_INFO_GENE
 * s -> struttura da copiare
 * 
 * Copia gli oggetti relativi alle informazione
 * sul gene.
 */
function copy_info_gene(s){
    
    var copy_info;
    copy_reg = {
    	"sequence_id" : s.sequence_id,
        "program_version" : s.program_version,
        "file_format_version" : s.file_format_version,
        "gene" : s.gene    
    };
    
    return copy_reg;    
}

/* INIT
 * 
 * Inizializza tutte le funzione per disegnare la struttura. 
 * Carica i dati dal file json relativo al gene selezionato.
 * Di default carica "ATP6AP1.json"
 */
function init(){
    
    //stringa per il pathname del file json
    var string = "json/";
    string = string.concat(s.property("value"), ".json");
    
    //console.log(string);
    //carica i dati contenuti nel file json e richiama le funzioni per disegnare la struttura
    //dell'isoforma
    d3.json(string, function(error, atp) {
	
	   console.log(error);
	   isoform = atp[0];
	   
	   original_info = copy_info_gene(isoform);
	   
	   //copia dell'array originale delle regioni
	   original_regions = copy_regions(isoform.regions);
	   //regioni
	   x = isoform_range(isoform.regions);
	   regions = regions_scaled(isoform.regions);
	
	   display_gene();
	
	   //boundaries
	   boundaries = isoform.boundaries;

	   //esoni
	   exons = isoform.exons;

	   //introni
	   introns = isoform.introns;
	   //console.log(introns);
	
	   //esoni, introni e boundaries ricostruiti
	   exons_restruct = exons_structure(exons, regions, boundaries);
	   introns_restruct = introns_structure(introns, regions, boundaries);
	   s_s_restruct = splice_site_structure(boundaries, regions);
	  
	   //disegna la struttura
	   line_i = draw_introns(svg_box, introns_restruct, x);
	   rect = draw_exons(svg_box, exons_restruct, x);
	   s_s = draw_splice_sites(svg_box, s_s_restruct, x);	  	
    });
}	