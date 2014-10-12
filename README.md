# Realizzazione di un’interfaccia grafica per la visualizzazione della struttura di un gene #

## Contesto biologico ##
Lo splicing alternativo è il meccanismo biologico che espande la complessità del proteoma 
negli organismi multicellulari. In virtù dello splicing alternativo, 
un gene può esprimere una molteplicità di trascritti (e quindi una molteciplità di proteine) 
combinando in maniera diversa i suoi esoni.

![fig_3.png](https://bitbucket.org/repo/5AzB7M/images/2408341880-fig_3.png)

L'interfaccia grafica è stata realizzata utilizzando la [libreria D3], al fine di fornire una visualizzazione della struttura in esoni e introni di un gene dovuta a splicing alternativo. 
[libreria D3]: http://d3js.org

## Requisiti ##
* Web server [Apache]
[Apache]: http://httpd.apache.org
* Browser con esecuzione di script Javascript attivata e compatibile con CSS BootStrap (Firefox consigliato)

## Funzionalità ##
L'interfaccia visualizza la struttura di un gene completa di esoni, introni e siti di splicing

![Schermata 2014-10-01 alle 11.32.03.png](https://bitbucket.org/repo/5AzB7M/images/2030944623-Schermata%202014-10-01%20alle%2011.32.03.png)

### Interattività ###
* Selezione ed espansione di elementi appartenenti ad una determinata regione
* Zoom della struttura genica
* Selezione e caricamento di altre strutture dati

## Upload delle strutture dati ##
Nella cartella **Json_file** sono presenti alcuni file JSON di esempio contenenti le strutture di un gene. Per caricare ulteriori strutture inserire i file nella cartella **Json_file** e il pathname nel file **config.json**.

## Note ##
La funzione di Zoom non è attiva se si utilizzano i browser Safari e Chrome a causa di un bug nell'esecuzione della funzione: 
```
#!javascript
d3.behavior.zoom()
```
che non permette lo zoom all'interno di una finestra creata con la classe SVG.
L'interfaccia è raggiungibile all'indirizzo: http://www.pen.statistica.unimib.it/~dellavedova/pinw