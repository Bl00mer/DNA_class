'''
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
    das einige Berechnung mit hilfe von Modelen durchgeführt werden, daher kann es sein das die Werte nicht richtig sind. 
    Außerdem kann dieses Programm nicht Spleißen, also kann nur Prokaryonten DNA oder Eukaryonten cDNA verwendet werden.
    Für die Suche nach Open Reading Frames werden werden 'ATG' als Start Codon und 'TAG','TAA' und 'TGA' als Stopcodons verwendet. Die Minimale Länge eines ORF liebt bei 300 Nukleotiden.
    Um die Ergebnisse plotten zu können muss matplotlib.pyplot als plt importiert werden
    An einer Stelle der Plot Methode der DNA_Sequenz Klasse kann random importiert werden, um zufällige Eigenschaften gegeieinander zu plotten


IDEEN für die Zukunft:
    isoelektronischer punkt
    stabilität 
    flexibilität
    sekundäre struktur merkmale
    Die Wahrscheinlichkeit einer funktionellen Veränderung durch die Analyse von SNPs (single-nucleotide polymorphisms) in der Sequenz
    Nucleotide Base Codes (IUPAC). DNA Sequenzen mit IUPAC-Nomenklatur zb. 'N' oder 'M' verwerten
    
Datum: 13.6.2023
Autor: H.M. Blum
'''

class DNA_Sequenz():
    
    '''
    Das ist die Klasse DNA_Sequenz. Zum erstellen eines Objektes muss eine DNA Sequenz und ein Name übergeben werden. 
    Sobald ein Objekt erstellt wird werden allerhand berechnungen durchgeführt um Eigenschaften der DNA Sequenz zu bestimmen. 
    Des weiteren gibt es Methoden (Show_distribution und Plot) die zur visualisierung verwendet werden können.
    
    
    
    Berechnete DNA Eigenschaften:
        revers_DNA: gibt den komplementären strang aus
        Länge: Länge der Nukleotidsequenz
        Base_Count: jaweilige anzahl der Basen A, T, G und C
        GC: GC gehalt
        Gewicht: Gewicht in Dalton
        copies_per_microG: Anzahl der Kopien pro Mikrogramm
        ORF_Finder: Erstellt mithilfe vorgegebener Suchparameter Objekte der ORF Klasse
        Codon_Usage: Errechnet den Codon Usage einer Sequnz und gibt eine Liste mit Anteilen aus
    Sonstige Funktionen:
        Show_Distribution: Zeigt die Trinucleotid Distribution als Heatmap
        Plot: damit kann man sich die Eigenschaften der ORFs eines DNA_Sequenz Klassenobjekts ausgeben lassen
    
    '''
    
    codons = [ # alle trinucleotid codons in alphabetisch reihenfolge
        'AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
        'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
        'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
        'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
        'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
        'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
        'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
        'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']
    
    instances = []# hier speichere ich alle erstellten Objekte der Klasse DNA_Sequenz
    
    def __init__(self, DNA, Name):
        self.Name = Name
        self.DNA = self.CHECK_UP(DNA) # Sequenz überprüfen
        self.revers_DNA = self.revers_DNA()
        
        self.Länge = len(self.DNA) 
        self.Base_Count=self.Base_Count()
        self.GC = sum(self.Base_Count[2:])/self.Länge # c und g sind die beiden letzten elemente in Base_Count--> summe g und c geteilt durch länge
        self.Gewicht=self.Gewicht()
        self.copies_per_microG = (6.022*10**23)/(self.Gewicht*10**6) # Avogadro[1/mol] /Seq-weight[g/mol]*10**6 = [Copies of Sequence/ng]
        
        DNA_Sequenz.instances.append(self)
        self.ORFs = self.ORF_Finder()
    def __str__(self):
        Info = (f'Name der DNA Sequenz: {self.Name},\n'+
                f'Gesamtlänge: {self.Länge}, '+
                f'GC:{self.GC*100:.2f}, '+
                f'Gewicht: {self.Gewicht} [g/mol], '+
                f'Copies per ug: {"{:.1e}".format(self.copies_per_microG)}, '+# damit formatiere ich das ergebnis in eine exponential schreibweise 
                f'Gefundene ORFs:{len(self.ORFs)}, '+
                f'ORF Dichte:{len(self.ORFs)/self.Länge}'
                "\n")
        return Info
    def __repr__(self):
        Info = (f'Name:{self.Name}'
                f' Länge in Basen: {self.Länge},'+
                f' Base in der Sequenz {self.Base_Count},'+
                f' GC:{self.GC},'+
                f' Gewicht [g/mol]: {self.Gewicht},'+
                f' Copies per µg {self.copies_per_microG},'+
                f' ORFs:{self.ORFs}, '+
                f'\n DNA:{self.DNA},'+
                f'\n revers DNA:{self.revers_DNA}'+
                '\n')

        return Info
    
    def CHECK_UP(self,DNA):
        '''
        Hier werden mögliche Zeilenumbrüche und Leerstellen in der Eingegebenen DNA Sequenz entfernt.
        Außerdem werden alle Buchstaben in Großbuchstaben umgewandelt
        
        '''
        DNA = DNA.replace('\n','')# entfernt die zeilenumbrüche
        DNA = DNA.replace(' ','')# entfernt von leerstellen
        DNA = DNA.upper()
        
        # """
        ### NBCI  ####   Nucleotide Base Codes (IUPAC) 
        
        # Mit der funktion will ich falls der eingegebene DNA code andere buchstaben (NBCI) enthält alternative dna sequencen erstellen. 
        
        # Das Problem ist nur die menge an sequencen die erstellt werden
        
        # Nucleotide Base Codes (IUPAC) (NBCI)
        #     Symbol: nucleotide(s)
        # """
        # base_dict = {
        #     'N': ['A', 'C', 'T', 'G'],
        #     'M': ['A', 'C'],
        #     'R': ['A', 'G'],
        #     'U': ['T'],
        #     'W': ['A', 'T'],
        #     'S': ['C', 'G'],
        #     'Y': ['C', 'T'],
        #     'K': ['G', 'T'],
        #     'V': ['A', 'C', 'G'],
        #     'H': ['A', 'C', 'T'],
        #     'D': ['A', 'G', 'T'],
        #     'B': ['C', 'G', 'T']}
        
        return DNA    
    
    def revers_DNA(self):
        '''
        Hier wird eine Komplementärer DNA strang in 5'-->3' Richtung gebildet.
        '''
        trans = self.DNA.maketrans("ATGC", "TACG") # umschreiben der dna von ATGC in TACG
        revers_DNA = self.DNA.translate(trans)[::-1] # übersezten und umdrehen der sequenz
        return revers_DNA
     
    def Base_Count(self):# erstellt eine liste mit den verwendeten 
        '''
        Zählt wie häufig jede Base in der Sequenz ist. Die Reihenfolge der ausgegebenen Liste ist A T G C 
        '''
        # Base_Count = [self.DNA.count('A') ,self.DNA.count('T'), self.DNA.count('G') ,self.DNA.count('C')] 
        
        Base_Count = [self.DNA.count('A') ,self.DNA.count('T'), self.DNA.count('G')] 
        Base_Count.append(self.Länge - sum(Base_Count))

        return Base_Count
    
    def Gewicht(self):
        
        '''
        Nucleic Acid Molecular Weight Conversions https://www.thermofisher.com/de/de/home/references/ambion-tech-support/rna-tools-and-calculators/dna-and-rna-molecular-weights-and-conversions.html
        Exact M.W. of ssRNA (e.g., RNA Transcript):
        M.W. = (An x 329.2) + (Un x 306.2) + (Cn x 305.2) + (Gn x 345.2) + 159ª
        An, Un, Cn, and Gn are the number of each respective nucleotide within the polynucleotide.
        M.W. calculated is valid at physiological pH.
        ªAddition of "159" to the M.W. takes into account the M.W. of a 5' triphosphate.
        
        Exact M.W. of ssDNA (e.g., Oligonucleotides):
        M.W. = (An x 313.2) + (Tn x 304.2) + (Cn x 289.2) + (Gn x 329.2) + 79.0ª
        An, Tn, Cn, and Gn are the number of each respective nucleotide within the polynucleotide.
        Addition of "79.0" to the M.W. takes into account the 5' monophosphate left by most restriction enzymes. No phosphate is present at the 5' end of strands made by primer extension.
        '''
        
        
        DNA_weight = (self.Base_Count[0] * 313.2) + (self.Base_Count[1] * 304.2) + (self.Base_Count[2] * 289.2) + (self.Base_Count[3] * 329.2) + 79.0
        return DNA_weight # g/mol
    
    def Codon_Usage(self,Sequenz):   #Trinucleotide distribution   
        '''
        Die Methode erstellt eine Liste mit dem Anteil des Codons an der gesamten Sequenz, in alphabetischer Reihenfolge.
        '''
        
        triplets = [Sequenz[i:i+3] for i in range(0, len(Sequenz), 3)]
        distribution = [0] * len(DNA_Sequenz.codons)  # Initialisiere die Verteilung mit Nullen
    
        for triplet in triplets:
            if triplet in DNA_Sequenz.codons:
                index = DNA_Sequenz.codons.index(triplet)
                distribution[index] += 1/(len(Sequenz)/3)
    
        return distribution
    
    def ORF_Finder(self, Startcodon=['ATG'], Stopcodon=['TAG','TAA','TGA'], Min_Länge = 300):
        '''
        Dieser Funktion wird eine DNA sequnec übergeben, aufgrund der gewählten suchparameter durchsucht 
        sie die Sequenc Frame für Frame nach nach einem Startcodon. Findet sie einen sucht sie nach einem 
        stopcodon im selben Frame. Wird zu dem Start ein passender Stopcodon gefunden wird eine 
        Instanz(also ein ORF bzw. Klassenobjekt) der Klasse ORF erstellt. Dann beginnt er ab dem letzten
        Stopcodon wieder die such nach einem neuen ORF, für jeden Frame und auch auf den reversen 
        komplementären Strang.
        
        die start und stop sequencen müssen derzeit 3 buchstaben lang sein
        
        
        '''

        ORF_Nummer = 0  
        ORFs=[]
        for o in [self.DNA,self.revers_DNA]:# zuerst wird der 5' -> 3' Strang und dann der reverse-komplementäre von 5' -> 3'
            for e in [0,1,2]:   # e ist der frame, oder das leseraster. 0 ist der Erste Frame
                i=0+e
                while i < self.Länge-2:####################         SOLLTE ICH HIER NICHT LEN DURCH die variable "Länge" ersetzen können--> weniger funktionsaufrufe??
                    if o[i:i+3] in Startcodon: 
                        j = i+3
                        while j < self.Länge-2:
                            if o[j:j+3] in Stopcodon:
                                    Frame = e+1    
                                    Länge_Seq = (j+3)-i
                                    DNA_Sequenz = o[i:j+3]
                                    if o==self.DNA:
                                        Strang = '+'
                                        Startposition = i
                                        Stopposition = j+3
                                    else:
                                        Strang = '-'
                                        Startposition = self.Länge-i-1
                                        Stopposition = self.Länge-j-2
                                    if Länge_Seq >= Min_Länge:
                                        ORF_Nummer +=1                            
                                        obj = ORF(ORF_Nummer,Strang,Frame, Startposition, Stopposition, Länge_Seq,DNA_Sequenz, self.Name)
                                        ORFs.append(obj) 
                                    i = j+3
                                    break
                            j += 3
                    i += 3           
        return ORFs         

    def Show_Distribution(self,Type=str):#,Type,CDS='' CDS sind ORFs von denen man weiß das sie codierend sind, und können sozusagen als key mit eingegeben werden
        '''
        Diese Funtkion visualisiert die Ergebnisse der Codon Usage von Objekten der DNA_Sequenz und der ORF Klasse.
        
        Wird 'ORF' als Type ausgewählt so wird eine liste mit den Codon Usages aller ORFs erstellt, die zu diesem DNA_Sequenz Klassenobjekt gehören.
        Wird 'DNA' als Type ausgewählt so wird eine liste mit den Codon Usage des DNA_Sequenz Klassenobjekts erstellt, 
        indem die DNA Sequenz in kleinere Sequenzen aufgeteilt und der Codon Usage dieser Sequenzen aufgelistet wird.

        Die Daten werden in einer HeatMap visualisiert, wobei ein Anteil von 1 100% darstellen würde. Daher geht die legende der heatmap nur bis 0,1 also 10%.
        '''
    
        label = f'{Type} {self.Name[:30]}'
        if Type == 'ORF':
            ergebnisse=[ORF.distribution for ORF in self.ORFs]# hiermit holt man die die eigenschaft distribution jedes ORF Objekts das in diesem Klassenobjekt gespeichert ist.
        elif Type == 'DNA':
            Schrittweite = 3000
            FensterGröße = 4500
            ergebnisse = [self.Codon_Usage(self.DNA[a-1:b]) for a, b in [(i, i+FensterGröße) for i in range(0, self.Länge - FensterGröße + 1, Schrittweite)]]

        else:
            print('Du musst schon angeben wovon du die distribuntion sehen willst!')
            return

        ergebnisse = list(zip(*ergebnisse))  # Transponieren der Ergebnismatrix

        # Erstellen des Heatmap-Diagramms
        # plt.imshow(ergebnisse, cmap='gist_rainbow', interpolation='nearest', aspect='auto')           # hier wird die legende automatisch erstellt, abhängig von den den codon usage werten
        # plt.imshow(ergebnisse, cmap='gist_rainbow', interpolation='nearest', aspect='auto', vmin=0, vmax=1) # hier geht die legende, von 0 bis 1. (unterschiede sind schwer zu erkennen)
        plt.imshow(ergebnisse, cmap='gist_rainbow', interpolation='nearest', aspect='auto', vmin=0, vmax=0.1) 
        
        
        # Title
        plt.title(label, fontdict=None, loc='center', pad=None)

        # Einstellen der Achsenbeschriftungen
        # plt.xticks(range(len(ergebnisse)), ['CDS {}'.format(i+1) for i in range(len(ergebnisse))], rotation=90)
        plt.xticks([])

        # plt.yticks(range(len(DNA_Sequenz.codons)), DNA_Sequenz.codons)
        plt.yticks([])

        # Hinzufügen von Farbleisten für die Farbskala
        plt.colorbar()

        # Anzeigen des Diagramms
        plt.show()
        
        return

    def Plot(self, Type=str, X=None, Y=None):
        '''
        Funktion zum Vergleichen von zwei Klasseneigenschaften von Objekten der DNA_Sequenz oder der ORF-Klasse.
        '''
        ############################################### Zufallsgenerator 
        Eigenschaften = ['länge','GC','tm','DNA_weight','Amount']
        if X is None or Y is None: # falls für x oder y keine eingabe gemacht wurde werden zufällige eigenschaften ausgewählt
            if X is None and Y is None:
                print('Beide Random')
                X,Y = random.choice(Eigenschaften), random.choice(Eigenschaften)# wählt einen zufälligen eintrag in einer liste
            elif X is None:
                print('X Random')
                X = random.choice(Eigenschaften)
            elif Y is None:
                print('Y Random')
                Y = random.choice(Eigenschaften)
        ###########################################
        
        label = f'ORFs {self.Name[:30]}'

        # [ORF.distribution for ORF in DNA_Sequenz.instances[0].ORFs]
        wertA_liste = [getattr(ORF, X) for ORF in self.ORFs]# The getattr() function returns the value of the specified attribute from the specified object.
        wertB_liste = [getattr(ORF, Y) for ORF in self.ORFs]

        # Title
        plt.title(label, fontdict=None, loc='center', pad=None)
        plt.scatter(x=wertA_liste, y=wertB_liste)
        plt.xlabel(X)
        plt.ylabel(Y)
        plt.show() 
                
class ORF():           #Um die Vererbung und die Erbreihenfolge festzulegen, muss nur die erbende Klasse im class-Aufruf in Klammern die allgemeine Klasse notiert werden:
    '''
    Ein objekt dieser Klasse ist ein Open Reading Frame (ORF). Zum erstellen eines Klassenobjekt braucht die 
    Klasse eine DNA Sequenz (str) und eine Startposition (int), da zur Bestimmung dieser variablen aber 
    auch Länge, Frame und Stopposition gehören werden sie gleich mit übergeben.
    Die Klasse hat einige Methoden welche die Nukleotid Sequenz des ORF nutzen um verschiedene Eigenschaften der DNA-Sequnez und der entsprechenden Aminosäure Sequenz zu berechnen.
    Berechnete DNA-Sequenz Eigenschaften:
        Base_Count*: Jeweilige Anzahl der basen A, T, G und C
        GC*: GC-Gehalt 
        Gewicht* **: das Gewicht der Nukleotid Sequenz
        copies_per_microG: Gibt an wieviele Kopien dieser DNA_Sequnez einen Mikrogramm wiegen
        tm**: die Schmelztemperatur des entsprechenden DNA-Doppelstrangs
    Berechnete Aminosäure-Sequenz Eigenschaften:
        AS: die Aminosäure Sequenz 
        AS_Länge: Die Länge der Aminosäure Sequenz
        Aroma**: die Aromatizität der Aminosäure Sequenz
        AS_weight: Das Gewicht der Aminosäure Sequenz in Dalton
        
    Zur Berechnung dieser Eigenschaft wurde eine Methode aus der DNA_Sequenz Klasse verwendet *
    Zur Berechnung dieser Eigenschaft wurde eine Mathematisches Model verwendet verwendet **
    '''
    instances = [] # die instanzen sind diese '<__main__.MyClass object at 0x0000029FF3F304C0>' kryptischen dinger. so kann ich alle klassen objekte die erstellt wurden in der klasse selbst speichern
    def __init__(self,ORF_Nr, Strang, Frame, Startposition, Stopposition, Länge, DNA, Name):# könnte ich hier nicht alles als alls liste übertragen? dann könnte ich die einzelnen eigenschaften über indicierung wieder rausholen
        self.ORF_Nr = ORF_Nr

        self.strang = Strang
        self.frame = Frame
        self.start = Startposition
        self.stop = Stopposition
        self.Länge = Länge
        self.sequence = DNA 
        self.DNA = self.sequence

        ORF.instances.append(self)
        
        ### DNA
        self.Base_Count= DNA_Sequenz.Base_Count(self)
        self.GC = sum(self.Base_Count[2:])/self.Länge    # hier benutze ich eine methode der Eltern klasse
        self.Gewicht = DNA_Sequenz.Gewicht(self)
        self.copies_per_microG = (6.022*10**23)/(self.Gewicht*10**6)
        
        self.tm = self.Schmelztemperatur()
        self.distribution=DNA_Sequenz.Codon_Usage(self, self.DNA)      
        # self.distribution=self.codon_usage()                   
             
        
        ### PROTEIN
        self.AS = self.translate()
        self.AS_Länge = len(self.AS)
        self.Aroma = self.calculate_aromaticity()
        self.AS_weight = self.calculate_protein_weight()
        # self.iP = self.isoelectronic_point()               

    def __str__(self):
        '''
        gibt eine String-Darstellung der Klasse zurück, die für den Benutzer 
        gedacht ist und zum Anzeigen von Informationen auf der Benutzeroberfläche verwendet 
        werden kann.
        Während __str__ für die Darstellung des Objekts in einer menschenlesbaren Form zuständig 
        ist, ist __repr__ für die Darstellung des Objekts in einer Form zuständig, die in 
        Python-Code verwendet werden kann.
        '''
     
        Info = (f'ORF Nr.{self.ORF_Nr:<3},'+
                f' Länge:{self.Länge:<5},'+
                f' Start:{self.start+1:<5},' + # wenn ich hier ',' benutze würde ich eine tulpe mit mehreren elementen als 'info' erhalten, so hängt alles zusammen
                f' Stop:{self.stop:<5},'+
                f' GC:{self.GC*100:.1f},'+
                f'DNA-Seq:{self.DNA[0:6]}...{self.DNA[-7:]} '+
                f' Codon Usage:{self.distribution},'+

                '\n')
        return Info
    
    def __repr__(self):
        Info = (
                f'ORF Nr.{self.ORF_Nr}\n'+
                f' Strang:{self.strang},'+
                f' Frame:{self.frame},'+
                f' Start:{self.start+1:<5},' + # wenn ich hier ',' benutze würde ich eine tulpe mit mehreren elementen als 'info' erhalten, so hängt alles zusammen
                f' Stop:{self.stop:<5},'+
                f' Länge:{self.Länge:<5},'+
                f' Basenanzahl [A,T,G,C]:{self.Base_Count},'+
                f' GC:{self.GC*100:.3f},'+
                f' DNA-Seq:{self.DNA},'+
                f' Tm[°C]:{round(self.tm,2):<6},'+
                f' Codon Usage:{self.distribution},'+
                f' Copies per ug: {"{:.1e}".format(self.copies_per_microG)},'+# damit formatiere ich das ergebnis in eine exponential schreibweise 
                f' AS Sequence:{self.AS},'+
                f' AS Länge:{self.AS_Länge},'+
                f' Aroma:{self.Aroma},'+
                f' AS Weight:{self.AS_weight},'+
                # f' Iso.ele.Point:{self.iP},'+
                '\n')
        return Info
    
    ### DNA #####################################################################
  
    def Schmelztemperatur(self):
        '''Wie wird die Schmelztemperatur TM eines Oligos berechnet?
        Die Schmelztemperatur TM charakterisiert die Stabilität eines DNA-Hybrids, das von einem Oligo 
        und seinem komplementären Strang gebildet wird. Bei TM 50% bindet ein Oligo an seinem komplementären Strang.
        
        Nutzen Sie unser Oligo Analyse Tool zur Berechnung der Schmelztemperatur. Ansonsten können folgende Formeln
        verwendet werden:
        
        Sequenz < 15 Basen
        TM [°C] = 2(nA + nT) + 4(nG + nC)
        
        Sequenz > 15 Basen
        TM [°C] = 69.3 + [41(nG + nC) / s – (650 / s)]
        
        n = Anzahl des Nukleosids X
        s = Menge aller Nukleoside der Sequenz
        
        Beispiel:
         Sequenz: GAA ATG AGT GCT CAT CAC TAC TTC CGC (27mer)
        nA=7; nC=8; nG=5; nT=7 s= 27
        
        Sequenz < 15 Basen
        TM [°C] = 2(7 + 7) + 4(5 + 8) = 80.0
        
        Sequenz > 15 Basen
        TM [°C] = 69.3 + [41(5 + 8) / 27] – (650 / 27) = 69.3 + 19.7 – 24.0 = 65.0
        
        
        es gibt noch alternative 
        
        Nächste-Nachbarn-Methode (https://www.sigmaaldrich.com/DE/de/technical-documents/protocol/genomics/pcr/oligos-melting-temp)
            '''
        tm  = 69.3+(41*(self.DNA.count('G')+self.DNA.count('C'))/self.Länge-(650/self.Länge)) # für sequnzen zwischen 5 und 200 bp(laut internetquelle)
        return tm #[°C]

    ### PROTEIN ################################################################ 
    def translate(self):
        codon_table = { #Es gibt insgesamt 64 (4^3) mögliche Triplet-Konformationen von den Basen A, T, G und C
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
        # '''Hier wird zuerst eine leere Liste protein erstellt. In der Schleife werden die Aminosäuren nacheinander zur Liste hinzugefügt. Am Ende wird die 
        # Liste mithilfe von join zu einem String zusammengefügt.'''
        AS = []
        for i in range(0, self.Länge, 3):
            codon = self.DNA[i:i+3]
            if codon in codon_table:
                amino_acid = codon_table[codon]
                AS.append(amino_acid)
        
        AS = ''.join(AS)
        return AS  

    def calculate_aromaticity(self):
        """
        Laut biopython muss nur die relative häufigkeit von 
        Phe+Thr+Tyr bestimmt werden
        Arome = Anzahl von F,W und Y geteilt dur die länge der AS-Seq
        
        """
        Aroma = (self.AS.count('F') + self.AS.count('W') + self.AS.count('Y'))/self.AS_Länge
        return Aroma
    
    def calculate_protein_weight(self): 
    # Massentabelle der Aminosäuren in Dalton
        aa_masses = {'A': 71.07, 'R': 156.18, 'N': 114.08, 'D': 115.08, 'C': 103.14,
                      'E': 129.11, 'Q': 128.13, 'G': 57.05, 'H': 137.14, 'I': 113.16,
                      'L': 113.16, 'K': 128.18, 'M': 131.20, 'F': 147.18, 'P': 97.12,
                      'S': 87.08, 'T': 101.11, 'W': 186.22, 'Y': 163.18, 'V': 99.13}
        AS_weight = sum(aa_masses.get(aa, 0) for aa in self.AS)
        return AS_weight        

    def isoelectronic_point(self): ########  NOT DONE
        '''
        beispiel
        Das berechnen des pI eines Proteins ist ein bisschen tricky, aber nicht kompliziert.
        Im grunde schreibt man das protein in etwa so auf:
        H3N-Asp-Gly-Glu-COO
        Dann schreibt man zu jedem die pKa Werte auf:
        9,60-3,65-/-4,25-2,19
        Dann beginnt man beliebige pH werte auzuprobieren und schaut wie sich die Ladungen bei dem gewählten pH ausgleichen.
        Bei pH7:(+)-(-)-(/)-(-)-(-)→ 1x(+)+3x(-)=2x(-)
        Bei pH10:(/)-(-)-(/)-(-)-(-) → 0x(+)+3x(-)=3x(-)
        Bei pH4:(+)-(-)-(/)-(/)-(-) → 1x(+)+2x(-)=1x(-)
        Bei pH3:(+)-(/)-(/)-(/)-(-) → 1x(+)+1x(-)=0
        Bei einem pH Wert von 3 ist die Nettoladung gleich 0. Um den pI zu erhalten nimmt man jetzt die pKa die als nächstes über und unter 3 liegen.
        (pKa<3+pKa>3)/2 =(2,19+3,65)/2 = pI= 2,92
        '''
        amino_acids = { # pka1 pka2 
            'A': (2.35, 9.87, None, 0),
            'R': (2.18, 9.04, 12.48, 0),####
            'N': (2.18, 8.72, None, 0),
            'D': (1.88, 9.60, 3.65, 0),####
            'C': (1.71, 10.78, 8.18, 0),###
            'Q': (2.17, 9.13, None, 0),
            'E': (2.19, 9.67, 4.25, 0),####
            'G': (2.34, 9.60, None, 0),
            'H': (1.78, 9.17, 6.00, 0),#####
            'I': (2.32, 9.76, None, 0),
            'L': (2.33, 9.74, None, 0),
            'K': (2.16, 9.06, 10.53, 0),######
            'M': (2.28, 9.21, None, 0),
            'F': (1.83, 9.13, None, 0),
            'P': (1.99, 10.6, None, 0),
            'S': (2.21, 9.15, None, 0),
            'T': (2.09, 9.10, None, 0),
            'W': (2.38, 9.39, None, 0),
            'Y': (2.20, 9.11, 10.07, 0),#####
            'V': (2.39, 9.74, None, 0)}

        pkas = []
        pkas.append(amino_acids[self.AS[0]][0])
        
        for i in self.AS[:-1]:
            if amino_acids[i][2] != None :
                #print(amino_acids[i][2])
                pkas.append(amino_acids[i][2])
            
        pkas.append(amino_acids[self.AS[-2]][1])

        # pH = 
        # while True:
            # a = 5
            # pH+= 0.01
            # if a==True:
            #     break
        # for i in range(0,14):
            
        #     break
        
        #iso_ele_point = 3
        return pkas#iso_ele_point
  

def app():
    '''
    Das ist die Funktion die startet sobald das Programm direkt (__name__  ==  "__main__") gestartet wird, würde man das programm über 'import ORF_class' als Modul
    importieren (also indirekt bzw. __name__  !=  "__main__") so würde die Funktion nicht automatisch starten.
    Die Funktion fordert den User dazu auf eine DNA Sequenz einzugeben (wie den Inhalt einer Genome FASTA Datei von NCBI MIT '>Namen blablabla' Zeile!). 
    Gibt man 'Test' ein so kann ein Ordner aus dem selben Verzeichnis mit mehrere Datein (FASTA, FAS, text ect.) eingelesen werden. Der Ordner muss nur 'Testdatein' genannt werden.

    '''

    ############################### Hier wird die Eingabe bestimmt
    # EINGABE = "CBDB1"
    # EINGABE = "Test"
    EINGABE = input("Hier Sequenz eingeben:")
    if EINGABE == "Test": ############ Wird 'Test' eingegeben werden die Text datein aus dem Ordner 'Testdatein' eingelesen
        ### TEST SEQUENZEN LADEN
        print('Lade und bearbeite Testsequenzen...')
        Arbeits_Verzeichniss = os.getcwd()# original verzeichniss merken# zeigt aktuelle datenpfand
        # print(os.listdir(os.getcwd())) #zeigt inhalt des aktuellen verzeichnisses
        os.chdir('Testdatein') # jetzt wechsele ich in das verzeichnis des ordner 'Testdatein'
        # print(os.listdir(os.getcwd()))  
        for i, item in enumerate(os.listdir(os.getcwd())[:]): # die zahl bestimmt bis zum wievielten element die elemente des verzeichnisses eingelesen werden sollen
            item_path = os.path.join(os.getcwd(), item)
            # print("Datei:", item, item_path)
            with open(item_path, "r") as file:
                content = file.read()
                # Name = content.split('\n', 1)[0] # Inhalt bis zum ersten '\n'
                Name = content.split('\n', 1)[0][1:] # Inhalt bis zum ersten '\n' ### das '>' Zeichen zu big jeder zeile wird entfernt
                content = content.split('\n', 1)[1] 
                Objekt = DNA_Sequenz(content, Name)            
                print('.')
                # print(f'{i+1}. Objekt \n'+ f'{Objekt}')

        os.chdir(Arbeits_Verzeichniss) #ins orginalverzeichniss zurückkehren
        
    elif EINGABE == "CBDB1":
        Arbeits_Verzeichniss = os.getcwd()
        
        os.chdir('Testdatein') 
        with open('Dehalococcoides mccartyi CBDB1.fasta', "r") as file:
            content = file.read()
            # Name = content.split('\n', 1)[0] # Inhalt bis zum ersten '\n'
            Name = content.split('\n', 1)[0][1:] # Inhalt bis zum ersten '\n' ### das '>' Zeichen zu beginn jedes eintrags wird entfernt
            content  = content.split('\n', 1)[1] 
            CBDB1 = DNA_Sequenz(content, Name)
        os.chdir(Arbeits_Verzeichniss) #ins orginalverzeichniss zurückkehren
    else:         
        Name = EINGABE.split('\n', 1)[0][1:] # Inhalt bis zum ersten '\n' ### das '>' Zeichen zu big jeder zeile wird entfernt
        content  = EINGABE.split('\n', 1)[1]
        DNA_Sequenz(content, Name)

                
    ######## HIER KANN MAN SICH ZB. VERSCHIEDENE WERTE AUSGEBEN LASSEN
    for i in range (len(DNA_Sequenz.instances)): 
        print(f'{i+1}. Objekt')
        print(str(DNA_Sequenz.instances[i])) # Hiermit kann man sich eine benutzerfreundliche Ansicht des KlassenObjektes (DNA_Sequenz) ausgeben lassen
        # print(repr(DNA_Sequenz.instances[i])) # Hiermit kann man sich die gesamten Informatione des Klassenobjekts ausgeben lassen

    print(
        # dir(DNA_Sequenz),
        # dir(ORF),
        # DNA_Sequenz.instances[0].ORFs,
        # len(ORF.instances),
        # DNA_Sequenz.instances[0].ORFs[0],
        # DNA_Sequenz.instances[0].ORFs[0].distribution,
        # len([ORF.distribution for instance in DNA_Sequenz.instances for ORF in instance.ORFs]),
        # [ORF.distribution for ORF in DNA_Sequenz.instances[0].ORFs],
        '')        
    # help(ORF)# HIFLE???
    # help(DNA_Sequenz)
    ####################################### HIER WIRD GEPLOTTED  ### Die Plot Plotfunktionen können nur im Python interpreter ausgeführt werden (nicht über die windows Konsole)
    try:
        for i in range (len(DNA_Sequenz.instances)):
            DNA_Sequenz.instances[i].Show_Distribution('DNA')
            DNA_Sequenz.instances[i].Show_Distribution('ORF')
            DNA_Sequenz.instances[i].Plot('ORF','Länge','GC')# Länge und GC sind hier nur Platzhalter und können durch andere Eigenschaften der ORF Klasse ersetzt werden
    except NameError:
        print('Die Plotfunktion ist in der Windoeskonsole nicht möglich.')
    
 
    

if __name__  ==  "__main__":
    import timeit
    import random
    import sys
    import os
    # je nachdenm ob das Programm über den interpreter oder die Windows konsole aufgerufen wurden, hat man verschiedene funktionen
    if sys.stdin.isatty():
        # Code, der ausgeführt wird, wenn der Code über die Windows-Konsole aufgerufen wird
        print("Code wird über die Windows-Konsole aufgerufen.")
        app()    
        for i in range (len(DNA_Sequenz.instances)): print(repr(DNA_Sequenz.instances[i]))# in der windowskonsole nochmal alle daten ausgeben

    else:
        # Code, der ausgeführt wird, wenn der Code im Interpreter ausgeführt wird
        import matplotlib.pyplot as plt
        print("Code wird im Interpreter ausgeführt.")
        app()    
    ### LAUFZEIT BESTIMMUNG  #########################################################
    # time = timeit.Timer(lambda:app())
    # print("Laufzeit der Anwendung:", time.timeit(number=10), "ms")
    
    
    # app()    
    ######## Ende der Anwendung    
    print('***')   