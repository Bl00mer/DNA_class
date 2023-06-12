# ORF_class
DNA und ORF Klasse 

Dieser Python-Code enthält zwei Klassen und eine Funktion.

Die erste Klasse heißt "DNA_Sequenz". Dieser Klasse wird eine DNA-Sequenz und ein Name übergeben. 
Anhand der Sequenz berechnet sie verschiedene Werte wie den GC-Gehalt, Codon Usage, Open Reading Frames usw.

Die zweite Klasse heißt "ORF". Objekte dieser Klasse werden automatisch erstellt, sobald ein Objekt der "DNA_Sequenz"-Klasse erstellt wird. 
Die Klasse berechnet einige Daten wie Schmelztemperatur, Gewicht usw. sowie Eigenschaften der Proteine, die sie möglicherweise kodieren könnten.

Die Funktion "app" wird gestartet, sobald dieser Code direkt ausgeführt wird, beispielsweise im Interpreter oder über die Windows-Konsole. 
Wenn dieses Programm als Modul importiert wird, wird die Funktion nicht automatisch ausgeführt. Diese Funktion fragt den Nutzer nach einer Eingabe, 
um ein Objekt der Klasse "DNA_Sequenz" zu erstellen.

Für genauere Erläuterungen siehe die Kommentare und Docstrings der entsprechenden Klassen und Funktionen.




Um dieses Programm (direkt) über die Windows Konsole auszuführen, musst du in das verzeichniss wechseln in dem du diesen code gespeichert hast (nutze dafür zb. den cd Befehl).
Dann gib die Folgende Zeile ein und bestätige:
    py -i ORF_class.py


BEACHTE:
    das einige Berechnung mit hilfe von Modelen durchgeführt werden, daher kann es sein das die Werte nicht immer richtig sind. 
    Außerdem kann dieses Programm nicht Spleißen, also kann nur Prokaryonten DNA oder Eukaryonten cDNA verwendet werden.
    Für die Suche nach Open Reading Frames werden werden 'ATG' als Start Codon und 'TAG','TAA' und 'TGA' als Stopcodons verwendet. Die Minimale Länge eines ORF liebt bei 300 Nukleotiden.
    Um die Ergebnisse plotten zu können muss matplotlib.pyplot als plt importiert werden.
    An einer Stelle der Plot Methode der DNA_Sequenz Klasse kann random importiert werden, um zufällige eigenschaften gegeieinander zu plotten.


IDEEN für die Zukunft:
    isoelektronischer punkt
    stabilität 
    flexibilität
    sekundäre struktur merkmale
    Die Wahrscheinlichkeit einer funktionellen Veränderung durch die Analyse von SNPs (single-nucleotide polymorphisms) in der Sequenz
    Nucleotide Base Codes (IUPAC). DNA Sequenzen mit IUPAC-Nomenklatur zb. 'N' oder 'M' verwerten
    
Datum: 13.6.2023
Autor: H.M. Blum
