import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;

public class Graph {
    private boolean adjMatrix[][];
    private int numVertices;

    
    public Graph(int numVertices) {
        this.numVertices = numVertices;
        adjMatrix = new boolean[numVertices][numVertices];
    }
    

    //Add edge
    public void addEdge(int i, int j) {
        adjMatrix[i][j] = true;
        adjMatrix[j][i] = true;
    }
    

    //Remove edge
    public void removeEdge(int i, int j) {
        adjMatrix[i][j] = false;
        adjMatrix[j][i] = false;
    }
    
    
    //Count number of edges for given row
    public int countEdges(int i) {
    	int sum = 0;
    	for(int j = 0; j < numVertices; j++) {
    		if(adjMatrix[i][j]) {
    			sum++;
    		}
    	}
    	return sum;
    }

    
    //Get the number of vertices
    public int getNumVertices() {
		return numVertices;
	}
    

    //Check if edge exists
	public boolean doesEdgeExist(int i, int j) {
		if(adjMatrix[i][j] == true) {
			return true;
		} else {
			return false;
		}
	}
	

	// Print the matrix
    public String toString() {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < numVertices; i++) {
            s.append(i + ": ");
            for (boolean j : adjMatrix[i]) {
                s.append((j ? 1 : 0) + " ");
            }
            s.append("\n");
        }
        
        return s.toString();
    }
    

    //Find groups which aren't connected. Returns arraylists of vertices. Each arraylist is a group which isn't connected to any other group
    public ArrayList<ArrayList<Integer>> toSubGraph() {
        ArrayList<Set<Integer>> result = new ArrayList<>();
        ArrayList<ArrayList<Integer>> resultArray = new ArrayList<>();
        
        //Changing table representation of network to sets of connected vertices
        for(int i = 0; i < numVertices - 1; i++){
            Set<Integer> subGraph = new TreeSet<>();
            subGraph.add(i);
            for(int j = i + 1; j < numVertices; j++){
                if(adjMatrix[i][j]){
                    subGraph.add(j);
                }
            }
            result.add(subGraph);
        }
        
	    //Merging sets to obtain subgraphs in network which aren't connected to each other 
	    boolean merge = true;
	    int numberOfMerges;
	    while (merge) {
	        numberOfMerges = 0;
	        for(int i = 0; i < result.size() - 1; i++){
	            for(int j = i + 1; j < result.size(); j++){
	                Set<Integer> subGraphCopy = new TreeSet<>(result.get(i));
	                subGraphCopy.retainAll(result.get(j));
	                if(!subGraphCopy.isEmpty()){
	                    result.get(i).addAll(result.get(j));
	                    result.remove(j);
	                    numberOfMerges++;
	                    j--;
	                }
	            }
	        }
	        if(numberOfMerges == 0){
	            merge = false;
	        }
	    }
	    
	    //Rewrite sets to arraylist
	    for(int i = 0; i < result.size(); i++){
	        ArrayList<Integer> subGraph = new ArrayList<>();
	        subGraph.addAll(result.get(i));
	        resultArray.add(subGraph);
	    }
        
        return resultArray;
    }
}
