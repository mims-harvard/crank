#ifndef snap_crank_h
#define snap_crank_h
#include "Snap.h"

/////////////////////////////////////////////////
/// CRank graph generator
class TCRank {
public:
  static void RndConnectInsideCommunity(PUNGraph& Graph, const TIntV& CmtyV, const double& Prob, TRnd& Rnd);
  static void AddConfModel(PUNGraph& Graph, const TIntV& DegSeqV, const int& NAdd, TRnd& Rnd);
  static void ProbaConnectInsideCommunity(PUNGraph& Graph, const TIntV& CmtyV, const TVec<TFltV>& ProbaVV, TRnd& Rnd);
  static PUNGraph GenProbaModel(TVec<TIntV>& CmtyVV, const THash<TIntPr,TFlt>& CProbaH, const THash<TIntPr,TFlt>& EdgeProbaH, const THash<TIntTr,TFlt>& CEdgeProba, const PUNGraph& OrigGraph, TRnd& Rnd=TInt::Rnd);
  static PUNGraph GenProbaModelAlpha(TVec<TIntV>& CmtyVV, const THash<TIntPr,TFlt>& CProbaH, const THash<TIntPr,TFlt>& EdgeProbaH, const THash<TIntTr,TFlt>& CEdgeProba, const PUNGraph& OrigGraph, const double Alpha, TRnd& Rnd=TInt::Rnd);
  static PUNGraph GenAGM(PUNGraph& Graph, TVec<TIntV>& CmtyVV, TRnd& Rnd=TInt::Rnd);
  static PUNGraph GenAGM(TVec<TIntV>& CmtyVV, const TFltV& CProbV, TRnd& Rnd, const double PNoCom = -1.0);
  static PUNGraph GenAGM(TVec<TIntV>& CmtyVV, const double& DensityCoef, const double& ScaleCoef, TRnd& Rnd);
  static PUNGraph GenAGMAlpha(PUNGraph& Graph, TVec<TIntV >& CmtyVV, const double Alpha = 0.15, TRnd& Rnd=TInt::Rnd);
  static PUNGraph GenAGMAlpha(TVec<TIntV>& CmtyVV, const TFltV& CProbV, const TIntV& DegSegV, const int& Edges, const double& Alpha, TRnd& Rnd, const double PNoCom = -1.0);
  static void GenCmtyProbV(TFltV& CProbV, TVec<TIntV>& CmtyVV, const double& DensityCoef, const double& ScaleCoef);
  static void AggrRank(TFltV& CRankV, TFltV& LklV, TFltV& DnsV, TFltV& BndV, TFltV& AllV, const double Prior = 0.05, int NBin = -1, const int MaxIter = 20, const double CorrStop = 1, const TFltV& SideV = TFltV());
  static void Prioritize(TFltV& CRankV, PUNGraph& Graph, TVec<TIntV >& CmtyVV, TRnd& Rnd=TInt::Rnd, const double Alpha = 0.15, const double Prior = 0.05, int NBin = -1, const int MaxIter = 20, const double CorrStop = 1, const bool FitProba = true, const THash<TIntPr,TFlt>& CProbaH = THash<TIntPr,TFlt>(), const THash<TIntPr,TFlt>& EdgeProbaH = THash<TIntPr,TFlt>(), const THash<TIntTr,TFlt>& CEdgeProbaH = THash<TIntTr,TFlt>(), const TFltV& SideV = TFltV());
  static void GetBayesFactors(TFltV& BayesFactorV, const TFltV& RankV, const TIntV& BinnedDataV, const int& NBin, const double& Prior);
  static void SmoothTukey(TFltV& InV);
};

/////////////////////////////////////////////////
/// CRank prioritization model
class TCRankMetric {
public:
  static void GetLikelihoodCmty(TFltV& LklV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV);
  static void GetDensityCmty(TFltV& DnsV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV);
  static void GetBoundaryCmty(TFltV& BndV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV);
  static void GetAllegianceCmty(TFltV& AllV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV);
  static void GetPriorMetr(TFltV& PmV, const TFltV& OrigV, const TFltV& PertV);
};

/////////////////////////////////////////////////
/// CRank utilities
class TCRankUtil {
public:
  static void GenPLSeq(TIntV& SzSeq,const int& SeqLen, const double& Alpha, TRnd& Rnd, const int& Min, const int& Max);
  static void GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const TVec<TIntV>& CmtyVV);
  static void GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const TVec<TIntV>& CmtyVV, const TIntV& NIDV);
  static void GetNodeMembership(TIntH& NIDComVH, const THash<TInt,TIntV >& CmtyVH);
  static void GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const TVec<TIntSet>& CmtyVV);
  static void GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const THash<TInt,TIntV >& CmtyVH);
  static void GetNodeMembership(THash<TInt,TIntV >& NIDComVH, const THash<TInt,TIntV >& CmtyVH);
  static void GetNodeMembership(THash<TInt,TIntV >& NIDComVH, const TVec<TIntV >& CmtyVV);
  static void LoadCmtyVV(const TStr& InFNm, TVec<TIntV>& CmtyVV);
  static void LoadCmtyVV(const TStr& InFNm, TVec<TIntV>& CmtyVV, TStrHash<TInt>& StrToNIdH, const int BeginCol, const int MinSz = 3, const TSsFmt Sep = ssfTabSep);
  static void DumpCmtyVV(const TStr& OutFNm, const TVec<TIntV>& CmtyVV);
  static void DumpCmtyVV(const TStr OutFNm, TVec<TIntV>& CmtyVV, TIntStrH& NIDNmH);
  static int TotalMemberships(const TVec<TIntV>& CmtyVV);
  static int Intersection(const TIntV& C1, const TIntV& C2);
  static void GetIntersection(const THashSet<TInt>& A, const THashSet<TInt>& B, THashSet<TInt>& C);
  static int Intersection(const THashSet<TInt>& A, const THashSet<TInt>& B);
  // Computes Modularity score of a set of nodes NIdV in a graph Graph
  static double GetModularity(const PUNGraph& Graph, const TIntV& NIdV, const int Edges);
  // Computes Modularity scores of a set of communities (each community is defined by its member nodes) in a graph Graph
  static void GetModularityCmty(TFltV& ModV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV);
  // Computes Conductance score of a set of nodes NIdV in a graph Graph
  static double GetConductance(const PUNGraph& Graph, const TIntSet& CmtyS, const int Edges);
  // Computes Conductance scores of a set of communities (each community is defined by its member nodes) in a graph Graph
  static void GetConductanceCmty(TFltV& PhiV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV);
  // Computes Random scores of a set of communities (each community is defined by its member nodes) in a graph Graph
  static void GetRandomCmty(TFltV& RndV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV, TRnd& Rnd=TInt::Rnd);
  // Assign ranks to data, dealing with ties appropriately, ranks begin at 1
  static void RankData(TFltV& OutV, const TFltV InV);
  static void CumMax(TFltV& DatV);
  static void ArgSort(TIntV& OutV, const TFltV InV);
  static double PearsonCorrelation(const TFltV& InA, const TFltV& InB);
  static void GetNbhCom(const PUNGraph& Graph, const int NID, TIntSet& NBCmtyS);
  static void SaveGephi(const TStr& OutFNm, const PUNGraph& G, const TVec<TIntV >& CmtyVVAtr, const double MaxSz, const double MinSz) {
    THash<TInt, TStr> TmpH;
    SaveGephi(OutFNm, G, CmtyVVAtr, MaxSz, MinSz, TmpH);
  }
  static void SaveGephi(const TStr& OutFNm, const PUNGraph& G, const TVec<TIntV >& CmtyVVAtr, const double MaxSz, const double MinSz, const THash<TInt, TStr>& NIDNameH) { 
    THash<TInt, TIntTr> TmpH; 
    SaveGephi(OutFNm, G, CmtyVVAtr, MaxSz, MinSz, NIDNameH, TmpH);
  }
  static void SaveGephi(const TStr& OutFNm, const PUNGraph& G, const TVec<TIntV >& CmtyVVAtr, const double MaxSz, const double MinSz, const THash<TInt, TStr>& NIDNameH, const THash<TInt, TIntTr >& NIDColorH);
  static void SaveBipartiteGephi(const TStr& OutFNm, const TIntV& NIDV, const TVec<TIntV>& CmtyVV, const double MaxSz, const double MinSz, const TIntStrH& NIDNameH, const THash<TInt, TIntTr >& NIDColorH, const THash<TInt, TIntTr >& CIDColorH);
  template <class PGraph>
  static PGraph LoadEdgeListStr(const TStr& InFNm, TIntStrH& NIDNameH, const int& SrcColId = 0, const int& DstColId = 1, const TSsFmt SsFmt = ssfTabSep) {
    TSsParser Ss(InFNm, SsFmt);
    PGraph Graph = PGraph::TObj::New();
    TStrHash<TInt> StrSet(Mega(1), true);
    while (Ss.Next()) {
      const int SrcNId = StrSet.AddKey(Ss[SrcColId]);
      const int DstNId = StrSet.AddKey(Ss[DstColId]);
      if (! Graph->IsNode(SrcNId)) { Graph->AddNode(SrcNId); }
      if (! Graph->IsNode(DstNId)) { Graph->AddNode(DstNId); }
      Graph->AddEdge(SrcNId, DstNId);
    }
    NIDNameH.Gen(StrSet.Len());
    for (int s = 0; s < StrSet.Len(); s++) { NIDNameH.AddDat(s, StrSet.GetKey(s)); }
    IAssert(NIDNameH.Len() == Graph->GetNodes());

    Graph->Defrag();
    return Graph;
  }
  template <class PGraph>
  static PGraph LoadEdgeListStr(const TStr& InFNm, TStrHash<TInt>& NodeNameH, const int& SrcColId = 0, const int& DstColId = 1, const TSsFmt SsFmt = ssfTabSep) {
    TSsParser Ss(InFNm, SsFmt);
    PGraph Graph = PGraph::TObj::New();
    TStrHash<TInt> StrSet(Mega(1), true);
    while (Ss.Next()) {
      const int SrcNId = StrSet.AddKey(Ss[SrcColId]);
      const int DstNId = StrSet.AddKey(Ss[DstColId]);
      if (! Graph->IsNode(SrcNId)) { Graph->AddNode(SrcNId); }
      if (! Graph->IsNode(DstNId)) { Graph->AddNode(DstNId); }
      Graph->AddEdge(SrcNId, DstNId);
    }
    NodeNameH = StrSet;
    NodeNameH.Pack();

    Graph->Defrag();
    return Graph;
  }

  template<class PGraph>
  static void GVizComGraph(const PGraph& Graph,const TVec<TIntV >& CmtyVV, const TStr& OutFNm, const TStr& Desc = TStr()){
    TStrV Colors = TStrV::GetV("red","blue","green","pink","cyan");
    TStrV Shapes = TStrV::GetV("ellipse","triangle","square","pentagon","hexagon");
    THash<TInt,TIntV> NIDComVH;
    GetNodeMembership(NIDComVH, CmtyVV);

    const TStr Ext = OutFNm.GetFExt();
    const TStr GraphFNm = OutFNm.GetSubStr(0, OutFNm.Len() - Ext.Len()) + "dot";
    const bool IsDir = HasGraphFlag(typename PGraph::TObj, gfDirected);
    FILE *F = fopen(GraphFNm.CStr(), "wt");
    if (! Desc.Empty()) fprintf(F, "/*****\n%s\n*****/\n\n", Desc.CStr());
    if (IsDir) { fprintf(F, "digraph G {\n"); } else { fprintf(F, "graph G {\n"); }
    fprintf(F, "  graph [splines=false overlap=false]\n"); //size=\"12,10\" ratio=fill
    fprintf(F, "  node  [width=0.3, height=0.3]\n");
    //nodes
    for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      int NID = NI.GetId();
      TIntV& CIDV = NIDComVH.GetDat(NID);
      IAssert(CIDV.Len() > 0);
      TStr ShapeNm = Shapes[(CIDV.Len()-1) % Shapes.Len()];
      TStr ColorNm = Colors[CIDV[0] % Colors.Len()];
      TStr NodeComLabel = TStr::Fmt("%d(",NID);
      for(int i=0;i<CIDV.Len();i++) {
        TStr TmpStr = TStr::Fmt("%d",int(CIDV[i]));
        NodeComLabel += TmpStr;
        if (i < CIDV.Len() - 1 ) { NodeComLabel += ","; }
      }
      NodeComLabel += ")";
      fprintf(F, "  %d [style=filled, shape=\"%s\" fillcolor=\"%s\" label=\"%s\"];\n", NI.GetId(), ShapeNm.CStr(), ColorNm.CStr(), NodeComLabel.CStr());
    }
  
    // edges
    for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
      if (NI.GetOutDeg()==0 && NI.GetInDeg()==0  ) { 
        fprintf(F, "%d;\n", NI.GetId()); }
      else {
        for (int e = 0; e < NI.GetOutDeg(); e++) {
          if (! IsDir && NI.GetId() > NI.GetOutNId(e)) { continue; }
          fprintf(F, "  %d %s %d;\n", NI.GetId(), IsDir?"->":"--", NI.GetOutNId(e)); 
        }
      }
    }
    if (! Desc.Empty()) {
      fprintf(F, "  label = \"\\n%s\\n\";", Desc.CStr());
      fprintf(F, "  fontsize=24;\n");
    }
    fprintf(F, "}\n");
    fclose(F);
    TSnap::TSnapDetail::GVizDoLayout(GraphFNm, OutFNm, gvlNeato);
  }
};

#endif
