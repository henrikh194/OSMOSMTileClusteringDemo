package parser;


/*
* Die Klasse OSMTileClusteringDemo kategorisiert einen beliebigen OSM-Datensatz
* in N Cluster ein die dann jeweils wiederum n-Quadranten der Länge nxn (meter) besitzen,
* z.B. 1 km x 1 km, außerdem wird jeder Cluster repräsentiert über Farbe und seine Top 3 Attribute.
* Die Ergebnisse werden im GeoJson-Format exportiert für die Verwendung in Visualisierungen,
* z.B. die Projektion auf die Ursprungskarte
*
*
*
* */




// Importieren der SAX XML Parsing Klassen
import org.xml.sax.*;
import org.xml.sax.helpers.DefaultHandler;

// Importieren für Parser, Dateioperationen, Collections und Streams
import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;
import java.io.File;
import java.io.FileWriter;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

// Import für JSON Export
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

public class OSMTileClusteringDemo {

    public static void main(String[] args) throws Exception {
        // Array von OSM-Dateien, die geparst werden sollen
        File[] osmFiles = {
                new File("C:/DUMMYPATH/HH.osm")
        };

        // Erstellen eines SAX-Handlers für OSM-Daten
        OSMHandler handler = new OSMHandler();

        // SAX Parser Factory erstellen
        SAXParserFactory factory = SAXParserFactory.newInstance();
        // SAX Parser erzeugen
        SAXParser parser = factory.newSAXParser();

        // Jede OSM-Datei parsen
        for (File f : osmFiles) {
            // Prüfen ob Datei existiert
            if (f.exists()) {
                // Ausgabe, welche Datei gerade geparst wird
                System.out.println("\nParsing file: " + f.getName());
                // Parsen der Datei mit dem Handler
                parser.parse(f, handler);
                // Ausgabe der bisher geparsten Nodes und Ways
                System.out.println("\rFinished parsing " + f.getName() + " (Nodes: " + handler.parsedNodes +
                        ", Ways: " + handler.parsedWays + ")");
            }
        }

        // Tiles erstellen, Größe 1 km
        handler.createTiles(1);

        // Anzahl der Cluster festlegen
        int k = 30;
        // Tiles in k Cluster gruppieren
        handler.clusterTiles(k);

        // Clusterinformationen ausgeben
        handler.printClusters();

        // Zufällige Tile abfragen und Clusterinformationen ausgeben
        handler.queryRandomTile();

        // GeoJSON Datei exportieren
        handler.exportTilesToGeoJSON("tiles_landuse.geojson");
    }

    static class OSMHandler extends DefaultHandler {

        // Map von Node-ID zu [Latitude, Longitude]
        private Map<Long, double[]> nodes = new HashMap<>();
        // Liste aller Ways
        private List<Way> ways = new ArrayList<>();
        // Temporäre Liste von Node-IDs für aktuellen Way
        private List<Long> currentWayNodes = new ArrayList<>();
        // Kategorie des aktuellen Ways
        private String currentCategory = null;
        // Flag, ob wir gerade einen Way parsen
        private boolean isWay = false;

        // Liste der erzeugten Tiles
        private List<Tile> tiles = new ArrayList<>();
        // Liste der Cluster
        private List<TileCluster> clusters = new ArrayList<>();

        // Fortschrittszähler für Nodes
        private long parsedNodes = 0;
        // Fortschrittszähler für Ways
        private long parsedWays = 0;

        /* Parsing */
        @Override
        public void startElement(String uri, String localName, String qName, Attributes attributes) {
            // Prüfen, ob das Element ein Node ist
            if ("node".equals(qName)) {
                // Node-ID auslesen
                long id = Long.parseLong(attributes.getValue("id"));
                // Latitude auslesen
                double lat = Double.parseDouble(attributes.getValue("lat"));
                // Longitude auslesen
                double lon = Double.parseDouble(attributes.getValue("lon"));
                // Node in Map speichern
                nodes.put(id, new double[]{lat, lon});

                // Zähler erhöhen
                parsedNodes++;
                // Fortschritt alle 100k Nodes ausgeben
                if (parsedNodes % 100_000 == 0) {
                    System.out.print("\rNodes parsed: " + parsedNodes + ", Ways parsed: " + parsedWays);
                }
            }

            // Prüfen, ob das Element ein Way ist
            if ("way".equals(qName)) {
                // Flag setzen
                isWay = true;
                // temporäre Node-Liste leeren
                currentWayNodes.clear();
                // Kategorie zurücksetzen
                currentCategory = null;

                // Zähler erhöhen
                parsedWays++;
                // Fortschritt alle 10k Ways ausgeben
                if (parsedWays % 10_000 == 0) {
                    System.out.print("\rNodes parsed: " + parsedNodes + ", Ways parsed: " + parsedWays);
                }
            }

            // Prüfen, ob Element ein nd innerhalb eines Ways ist
            if (isWay && "nd".equals(qName))
                // Node-ID des Ways speichern
                currentWayNodes.add(Long.parseLong(attributes.getValue("ref")));

            // Prüfen, ob Element ein tag innerhalb eines Ways ist
            if (isWay && "tag".equals(qName)) {
                // Schlüssel auslesen
                String k = attributes.getValue("k");
                // Wert auslesen
                String v = attributes.getValue("v");
                // Nur bestimmte Kategorien merken
                if ("landuse".equals(k) || "natural".equals(k) || "leisure".equals(k)) {
                    currentCategory = v;
                }
            }
        }

        @Override
        public void endElement(String uri, String localName, String qName) {
            // Prüfen, ob das Element ein Way beendet
            if ("way".equals(qName)) {
                // Wenn Kategorie vorhanden, Way speichern
                if (currentCategory != null) {
                    ways.add(new Way(new ArrayList<>(currentWayNodes), currentCategory));
                }
                // Flag zurücksetzen
                isWay = false;
            }
        }

        /* Tiles erstellen */
        public void createTiles(int kmSize) {
            // minimale Latitude ermitteln
            double minLat = nodes.values().stream().mapToDouble(p -> p[0]).min().getAsDouble();
            // maximale Latitude ermitteln
            double maxLat = nodes.values().stream().mapToDouble(p -> p[0]).max().getAsDouble();
            // minimale Longitude ermitteln
            double minLon = nodes.values().stream().mapToDouble(p -> p[1]).min().getAsDouble();
            // maximale Longitude ermitteln
            double maxLon = nodes.values().stream().mapToDouble(p -> p[1]).max().getAsDouble();

            // Schrittweite Latitude berechnen
            double latStep = kmSize * 0.00899;
            // Schrittweite Longitude berechnen
            double lonStep = kmSize * 0.015;

            // temporäre Map für Tiles
            Map<String, Tile> tileMap = new HashMap<>();
            // jeden Way durchgehen
            for (Way w : ways) {
                // Mittelpunkt des Ways berechnen
                double[] centroid = w.getCentroid(nodes);
                // Tile-Index Latitude berechnen
                int i = (int) ((centroid[0] - minLat) / latStep);
                // Tile-Index Longitude berechnen
                int j = (int) ((centroid[1] - minLon) / lonStep);
                // Schlüssel für Tile erzeugen
                String key = i + "_" + j;
                // Tile erstellen oder aus Map holen
                Tile t = tileMap.getOrDefault(key,
                        new Tile("Tile " + key, "DummyLocation " + key, i, j, minLat, minLon, latStep, lonStep));
                // Way zur Tile hinzufügen
                t.addWay(w, nodes);
                // Tile in Map speichern
                tileMap.put(key, t);
            }
            // alle Tiles in Liste speichern
            tiles.addAll(tileMap.values());
            // Anzahl der Tiles ausgeben
            System.out.println("\nTiles erstellt: " + tiles.size());
        }

        /* Clustering über K-Means */
        public void clusterTiles(int k) {
            // Feature-Vektor für jede Tile erstellen
            List<double[]> features = tiles.stream().map(Tile::getFeatureVector).collect(Collectors.toList());
            // K-Means Clustering durchführen

            int[] labels = KMeans.fit(features, k);
            double silhouette = KMeans.silhouetteScore(features, labels, k);
            System.out.println("Silhouette Score (k=" + k + "): " + String.format("%.4f", silhouette));

            // Cluster-Objekte erstellen
            Map<Integer, TileCluster> clusterMap = new HashMap<>();
            for (int i = 0; i < tiles.size(); i++) {
                // Label des Punktes
                int l = labels[i];
                // Cluster-Objekt erstellen falls noch nicht vorhanden
                clusterMap.putIfAbsent(l, new TileCluster(l));
                // Tile zum Cluster hinzufügen
                clusterMap.get(l).addTile(tiles.get(i));
            }
            // Cluster in Liste speichern
            clusters.addAll(clusterMap.values());

            // Top Features pro Cluster berechnen
            for (TileCluster c : clusters) c.computeTopFeatures();
        }

        public void printClusters() {
            // Ausgabe aller Cluster
            for (TileCluster c : clusters) System.out.println(c);
        }

        public void queryRandomTile() {
            // Prüfen, ob Tiles vorhanden
            if (tiles.isEmpty()) return;
            // Zufällige Tile auswählen
            Random rnd = new Random();
            Tile t = tiles.get(rnd.nextInt(tiles.size()));
            // Cluster-ID der Tile ermitteln
            int clusterId = clusters.stream().filter(c -> c.tiles.contains(t)).findFirst().map(c -> c.id).orElse(-1);
            System.out.println("\nRandom Tile: " + t.name + " / " + t.originLocation);
            // Cluster-Label ausgeben
            System.out.println("Wird zu Cluster " + clusterId + " zugeordnet: " +
                    (clusterId >= 0 ? clusters.stream().filter(c -> c.id == clusterId).findFirst().map(c -> c.label).orElse("N/A") : "N/A"));
        }

        /* GeoJSON Export */
        public void exportTilesToGeoJSON(String filename) throws Exception {
            // Liste für Features
            List<Map<String, Object>> features = new ArrayList<>();

            // Farben für Kategorien
            Map<String, String> categoryColors = new HashMap<>();
            categoryColors.put("water", "#1f78b4");
            categoryColors.put("lake", "#1f78b4");
            categoryColors.put("river", "#1f78b4");
            categoryColors.put("forest", "#33a02c");
            categoryColors.put("grass", "#b2df8a");
            categoryColors.put("park", "#33a02c");
            categoryColors.put("residential", "#ff7f00");
            categoryColors.put("industrial", "#e31a1c");
            categoryColors.put("commercial", "#6a3d9a");
            categoryColors.put("farmland", "#b2df8a");
            categoryColors.put("leisure", "#a6cee3");
            categoryColors.put("natural", "#b15928");
            categoryColors.put("meadow", "#c2f0c2");
            categoryColors.put("farm", "#fdbf6f");
            categoryColors.put("industrial_area", "#fb8072");

            // Jede Tile verarbeiten
            for (Tile t : tiles) {
                // Cluster-ID ermitteln
                int clusterId = clusters.stream().filter(c -> c.tiles.contains(t)).findFirst().map(c -> c.id).orElse(-1);
                // Cluster-Objekt abfragen
                TileCluster cluster = clusters.stream().filter(c -> c.id == clusterId).findFirst().orElse(null);

                // Feature Map erstellen
                Map<String, Object> feature = new HashMap<>();
                feature.put("type", "Feature");

                // Properties erstellen
                Map<String, Object> properties = new HashMap<>();
                properties.put("name", t.name);
                properties.put("originLocation", t.originLocation);
                properties.put("cluster", clusterId);
                properties.put("label", cluster != null ? cluster.label : "N/A");

                // Farbe basierend auf Hauptkategorie
                String color = "#808080";
                if (cluster != null) {
                    String mainCategory = cluster.getMainCategory();
                    if (mainCategory != null && categoryColors.containsKey(mainCategory)) {
                        color = categoryColors.get(mainCategory);
                    }
                }
                properties.put("color", color);

                feature.put("properties", properties);

                // Polygon-Koordinaten für Tile
                Map<String, Object> geometry = new HashMap<>();
                geometry.put("type", "Polygon");
                geometry.put("coordinates", List.of(List.of(
                        List.of(t.lon1, t.lat1),
                        List.of(t.lon2, t.lat1),
                        List.of(t.lon2, t.lat2),
                        List.of(t.lon1, t.lat2),
                        List.of(t.lon1, t.lat1)
                )));
                feature.put("geometry", geometry);

                // Feature zur Liste hinzufügen
                features.add(feature);
            }

            // FeatureCollection Map erstellen
            Map<String, Object> geojson = new HashMap<>();
            geojson.put("type", "FeatureCollection");
            geojson.put("features", features);

            // In Datei schreiben
            try (FileWriter writer = new FileWriter(filename)) {
                Gson gson = new GsonBuilder().setPrettyPrinting().create();
                gson.toJson(geojson, writer);
            }

            // Erfolgsmeldung
            System.out.println("GeoJSON exportiert: " + filename);
        }

        /* Hilfsklassen */
        static class Way {
            // Node-IDs des Ways
            List<Long> nodes;
            // Kategorie des Ways
            String category;

            Way(List<Long> nodes, String category) {
                this.nodes = nodes;
                this.category = category;
            }

            // Centroid berechnen
            double[] getCentroid(Map<Long, double[]> nodeMap) {
                double sumLat = 0, sumLon = 0;
                int count = 0;
                for (Long id : nodes) {
                    double[] p = nodeMap.get(id);
                    if (p != null) {
                        sumLat += p[0];
                        sumLon += p[1];
                        count++;
                    }
                }
                return new double[]{sumLat / count, sumLon / count};
            }

            // Fläche berechnen
            double computeArea(Map<Long, double[]> nodeMap) {
                List<double[]> points = new ArrayList<>();
                for (Long id : nodes) {
                    double[] p = nodeMap.get(id);
                    if (p != null) points.add(p);
                }
                if (points.size() < 3) return 0;
                double sum = 0;
                for (int i = 0; i < points.size(); i++) {
                    double[] p1 = points.get(i);
                    double[] p2 = points.get((i + 1) % points.size());
                    sum += (p1[1] * p2[0] - p2[1] * p1[0]);
                }
                return Math.abs(sum / 2.0);
            }
        }

        static class Tile {
            String name, originLocation;
            Map<String, Double> featureCounts = new HashMap<>();
            double lat1, lat2, lon1, lon2;
            static List<String> featureKeys = new ArrayList<>();

            Tile(String name, String originLocation, int i, int j, double minLat, double minLon, double latStep, double lonStep) {
                this.name = name;
                this.originLocation = originLocation;
                this.lat1 = minLat + i * latStep;
                this.lon1 = minLon + j * lonStep;
                this.lat2 = lat1 + latStep;
                this.lon2 = lon1 + lonStep;
            }

            // Way zur Tile hinzufügen
            void addWay(Way w, Map<Long, double[]> nodeMap) {
                if (w.category != null) {
                    double area = w.computeArea(nodeMap);
                    featureCounts.merge(w.category, area, Double::sum);
                    if (!featureKeys.contains(w.category)) featureKeys.add(w.category);
                }
            }

            // Feature-Vektor für K-Means
            double[] getFeatureVector() {
                double[] fv = new double[featureKeys.size()];
                for (int i = 0; i < featureKeys.size(); i++) {
                    fv[i] = featureCounts.getOrDefault(featureKeys.get(i), 0.0);
                }
                return fv;
            }
        }

        static class TileCluster {
            int id;
            List<Tile> tiles = new ArrayList<>();
            String label;

            TileCluster(int id) { this.id = id; }

            void addTile(Tile t) { tiles.add(t); }

            void computeTopFeatures() {
                Map<String, Double> sums = new HashMap<>();
                for (Tile t : tiles)
                    for (Map.Entry<String, Double> e : t.featureCounts.entrySet())
                        sums.merge(e.getKey(), e.getValue(), Double::sum);

                List<Map.Entry<String, Double>> top = sums.entrySet().stream()
                        .sorted(Map.Entry.<String, Double>comparingByValue().reversed())
                        .limit(3).collect(Collectors.toList());

                double total = sums.values().stream().mapToDouble(Double::doubleValue).sum();
                label = top.stream().map(e -> e.getKey() + " " + String.format("%.2f", e.getValue() / total))
                        .collect(Collectors.joining(", "));
            }

            String getMainCategory() {
                if (tiles.isEmpty()) return null;
                Map<String, Double> sums = new HashMap<>();
                for (Tile t : tiles)
                    for (Map.Entry<String, Double> e : t.featureCounts.entrySet())
                        sums.merge(e.getKey(), e.getValue(), Double::sum);
                return sums.entrySet().stream()
                        .max(Map.Entry.comparingByValue())
                        .map(Map.Entry::getKey)
                        .orElse(null);
            }

            public String toString() {
                return "Cluster " + id + ": " + label + " (" + tiles.size() + " Tiles)";
            }
        }
        // static Klasse für K-Means Clustering
            //iteratives "Zusammenstellen anhand der euklidischen Distanz
        static class KMeans {
            // Methode, um K-Means auf eine Liste von allen Feature-Vektoren anzuwenden
            // data: Liste von double-Arrays (Feature-Vektoren)
            // k: Anzahl der Cluster
            static int[] fit(List<double[]> data, int k) {
                // n = Anzahl der Datenpunkte (Tiles)
                int n = data.size();
                // dim = Dimension der Feature-Vektoren (Anzahl der Kategorien)
                int dim = data.get(0).length;
                // Zufallsobjekt für die Initialisierung der Zentroiden
                Random rnd = new Random();

                // Zentroiden-Array erstellen: k Cluster, jedes mit dim Dimensionen
                double[][] centroids = new double[k][dim];
                // Zentroiden zufällig aus den Daten initialisieren
                for (int i = 0; i < k; i++)
                    centroids[i] = data.get(rnd.nextInt(n)).clone(); // klonen, damit Originaldaten erhalten bleibt

                // Array für die Clusterzuweisungen jedes Datenpunkts
                int[] labels = new int[n];

                // Iterationen der K-Means Schleife (typisch 30 Iterationen, mehr möglich, aber Performance(!))
                for (int iter = 0; iter < 30; iter++) {
                    // jedem Punkt das nächste Zentroid zuordnen
                    for (int i = 0; i < n; i++) {
                        double minDist = Double.MAX_VALUE; // Startwert für die minimale Distanz
                        int best = 0; // Index des nächsten Zentroids
                        // Über alle Cluster-Zentroiden prüfen
                        for (int c = 0; c < k; c++) {
                            // Distanz zwischen Punkt und Zentroid berechnen
                            double d = dist(data.get(i), centroids[c]);
                            // Wenn Distanz kleiner als bisherige minimale Distanz
                            if (d < minDist) {
                                minDist = d;
                                best = c; // nächstes Cluster merken
                            }
                        }
                        // Clusterzuweisung speichern
                        labels[i] = best;
                    }

                    //Zentroiden neu berechnen (Mittelwert aller Punkte im Cluster)
                    double[][] newCentroids = new double[k][dim]; // Zwischenspeicher für neue Zentroiden
                    int[] counts = new int[k]; // Zähler, wie viele Punkte in jedem Cluster
                    for (int i = 0; i < n; i++) {
                        int l = labels[i]; // Cluster des Punktes
                        for (int d = 0; d < dim; d++)
                            newCentroids[l][d] += data.get(i)[d]; // Feature-Werte aufsummieren
                        counts[l]++; // Zähler erhöhen
                    }
                    // Durchschnitt berechnen, um neue Zentroiden zu setzen
                    for (int c = 0; c < k; c++)
                        if (counts[c] > 0) // Division nur, wenn Cluster Punkte enthält
                            for (int d = 0; d < dim; d++)
                                newCentroids[c][d] /= counts[c];
                    // Zentroiden auf die neuen Werte setzen
                    centroids = newCentroids;
                }
                // Clusterlabels für jeden Datenpunkt zurückgeben
                return labels;
            }

            // Hilfsmethode: Euklidische Distanz zwischen zwei Vektoren berechnen
            static double dist(double[] a, double[] b) {
                double sum = 0;
                for (int i = 0; i < a.length; i++)
                    sum += (a[i] - b[i]) * (a[i] - b[i]); // Differenz quadrieren und aufsummieren
                return Math.sqrt(sum); // Quadratwurzel ziehen = Euklidische Distanz
            }


            static double silhouetteScore(List<double[]> data, int[] labels, int k) {
                int n = data.size();

                // Cluster → Indizes der Punkte
                Map<Integer, List<Integer>> clusters = new HashMap<>();
                for (int i = 0; i < n; i++) {
                    clusters.computeIfAbsent(labels[i], x -> new ArrayList<>()).add(i);
                }

                double totalScore = 0.0;
                // Berechnung des Silhouettenkoeffizienten für alle Datenpunkte:
                //  Für jeden Punkt i wird die mittlere Intra-Cluster-Distanz a(i)
                // sowie   die minimale mittlere Inter-Cluster-Distanz b(i) bestimmt.
                // Aus   diesen beiden Werten wird der Shilhhouettenwert s(i) berechnet,
                // der angibt,  wie präzise der Punkt seinem zugewiesenen Cluster folgt.
                // Der  GesamteScore ergibt sich als Mittelwert  über ale s(i).
                for (int i = 0; i < n; i++) {
                    int clusterId = labels[i];
                    double[] point = data.get(i);

                    // a(i): mittlere Distanz zu Punkten im selben Cluster
                    double a = 0.0;
                    List<Integer> sameCluster = clusters.get(clusterId);

                    if (sameCluster.size() > 1) {
                        for (int j : sameCluster) {
                            if (j != i) {
                                a += dist(point, data.get(j));
                            }
                        }
                        a /= (sameCluster.size() - 1);
                    }

                    // b(i): minimale mittlere Distanz zu einem anderen Cluster
                    double b = Double.MAX_VALUE;

                    for (Map.Entry<Integer, List<Integer>> entry : clusters.entrySet()) {
                        int otherClusterId = entry.getKey();
                        if (otherClusterId == clusterId) continue;

                        double distSum = 0.0;
                        for (int j : entry.getValue()) {
                            distSum += dist(point, data.get(j));
                        }
                        double meanDist = distSum / entry.getValue().size();
                        b = Math.min(b, meanDist);
                    }

                    // Silhouette für Punkt i
                    double s;
                    if (a == 0 && b == 0) {
                        s = 0;
                    } else {
                        s = (b - a) / Math.max(a, b);
                    }

                    totalScore += s;
                }

                return totalScore / n;
            }

        }

    }
}
