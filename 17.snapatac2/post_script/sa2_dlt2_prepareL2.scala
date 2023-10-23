import scala.io.Source
import scala.collection.immutable._
import java.io._

case class Cell(barcode: String, clusterId: Int)
case class ClusterSum(head: Tuple2[String, String], cells: List[Cell])

def loadClusterSumFromFile(file: String, sep:String = ","): ClusterSum = {
  val lines= Source.fromFile(file).getLines.toList
  val firstline= lines.head.trim.split(sep)
  val head: Tuple2[String, String] = (firstline(0), firstline(1))
  val cells = lines.tail.map(x => x.trim.split(sep)).map(x => Cell(x(0), x(1).toInt))
  ClusterSum(head, cells)
}

def getClusterSize(x: ClusterSum): Map[Int, Int] = {
  x.cells.groupBy(_.clusterId).map(t => (t._1, t._2.length))
}

def writeMap2csv(fnm: String, lines: Map[Int,Int]): Unit = {
  val file = new File(fnm)
  val bw = new BufferedWriter(new FileWriter(file))
  lines.foreach(t => bw.write(s"${t._1},${t._2}\n"))
  bw.close()
}

@main def Hello(params: String*):Unit = {
  val clusterSum = loadClusterSumFromFile(
    file = "../result/clustering_sum_L1/sa2_L1_r0.4_barcodes2id.csv", sep = ",")
  val cluster2size = getClusterSize(clusterSum)
  writeMap2csv(fnm = "../resource/sa2_dlt2_L1_cluster2size.csv", lines = cluster2size)
}
