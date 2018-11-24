#include "stdafx.h"
#include "Snap.h"
#include "crank.h"

/////////////////////////////////////////////////
// CRank network model (uses Affiliation Graph Model (AGM) graph generator)

void TCRank::Prioritize(TFltV& CRankV, PUNGraph& Graph, TVec<TIntV >& CmtyVV, TRnd& Rnd, const double Alpha, const double Prior, int NBin, const int MaxIter, const double CorrStop, const bool FitProba, const THash<TIntPr,TFlt>& CProbaH, const THash<TIntPr,TFlt>& EdgeProbaH, const THash<TIntTr,TFlt>& CEdgeProbaH, const TFltV& SideV) {
  PUNGraph OrigG, PertG;
  if (FitProba) {
    OrigG = TCRank::GenAGM(Graph, CmtyVV, Rnd);
    PertG = TCRank::GenAGMAlpha(Graph, CmtyVV, Alpha, Rnd);
  }
  else {
    OrigG = TCRank::TCRank::GenProbaModel(CmtyVV, CProbaH, EdgeProbaH, CEdgeProbaH, Graph, Rnd);
    PertG = TCRank::GenProbaModelAlpha(CmtyVV, CProbaH, EdgeProbaH, CEdgeProbaH, Graph, Alpha, Rnd);
  }
  printf("G: %d G(alpha = 0): %d, G(alpha = %f): %d\n", Graph->GetEdges(), OrigG->GetEdges(), Alpha, PertG->GetEdges());
  // structural features in original network
  printf("Structural features in G(alpha = 0)\n");
  TFltV LklOrigV, DnsOrigV, BndOrigV, AllOrigV;
  TCRankMetric::GetLikelihoodCmty(LklOrigV, OrigG, CmtyVV);
  TCRankMetric::GetDensityCmty(DnsOrigV, OrigG, CmtyVV);
  TCRankMetric::GetBoundaryCmty(BndOrigV, OrigG, CmtyVV);
  TCRankMetric::GetAllegianceCmty(AllOrigV, OrigG, CmtyVV);
  // structural features in perturbed network
  printf("Structural features in G(alpha = %f)\n", Alpha);
  TFltV LklPertV, DnsPertV, BndPertV, AllPertV;
  TCRankMetric::GetLikelihoodCmty(LklPertV, PertG, CmtyVV);
  TCRankMetric::GetDensityCmty(DnsPertV, PertG, CmtyVV);
  TCRankMetric::GetBoundaryCmty(BndPertV, PertG, CmtyVV);
  TCRankMetric::GetAllegianceCmty(AllPertV, PertG, CmtyVV);
  // prioritization metrics
  printf("Prioritization metrics\n");
  TFltV LklV, DnsV, BndV, AllV;
  TCRankMetric::GetPriorMetr(LklV, LklOrigV, LklPertV);
  TCRankMetric::GetPriorMetr(DnsV, DnsOrigV, DnsPertV);
  TCRankMetric::GetPriorMetr(BndV, BndOrigV, BndPertV);
  TCRankMetric::GetPriorMetr(AllV, AllOrigV, AllPertV);
  // rank aggregation
  printf("Rank aggregation\n");
  TCRank::AggrRank(CRankV, LklV, DnsV, BndV, AllV, Prior, NBin, MaxIter, CorrStop, SideV);
}

void TCRank::AggrRank(TFltV& CRankV, TFltV& LklV, TFltV& DnsV, TFltV& BndV, TFltV& AllV, const double Prior, int NBin, const int MaxIter, const double CorrStop, const TFltV& SideV) {
  int Nr = LklV.Len();
  if (NBin == -1) {
    //NBin = fmax(3, (int) Nr / 50);
    NBin = 3;
    if (NBin < (int) Nr / 50) {
      NBin = (int) Nr / 50;
    }
  }
  printf("NBin: %d\n", NBin);
  //int Nrp = (int) fmax(1, Nr * Prior);
  int Nrp = 1;
  if (Nrp < Nr * Prior) {
    Nrp = (int) (Nr * Prior);
  }
  printf("Nrp: %d\n", Nrp);
  TFltV LklRankV; TIntV LklBinnedDataV;
  LklRankV.Gen(Nr); LklBinnedDataV.Gen(Nr);
  TCRankUtil::RankData(LklRankV, LklV);
  TFltV DnsRankV; TIntV DnsBinnedDataV;
  DnsRankV.Gen(Nr); DnsBinnedDataV.Gen(Nr);
  TCRankUtil::RankData(DnsRankV, DnsV);
  TFltV BndRankV; TIntV BndBinnedDataV;
  BndRankV.Gen(Nr); BndBinnedDataV.Gen(Nr);
  TCRankUtil::RankData(BndRankV, BndV);
  TFltV AllRankV; TIntV AllBinnedDataV;
  AllRankV.Gen(Nr); AllBinnedDataV.Gen(Nr);
  TCRankUtil::RankData(AllRankV, AllV);
  TFltV ExtraRankV; TIntV ExtraBinnedDataV;
  ExtraRankV.Gen(Nr); ExtraBinnedDataV.Gen(Nr);
  TCRankUtil::RankData(ExtraRankV, SideV);
  for (int i = 0; i < Nr; i++) {
    LklBinnedDataV[i] = (int) ceil(LklRankV[i] / Nr * NBin);
    DnsBinnedDataV[i] = (int) ceil(DnsRankV[i] / Nr * NBin);
    BndBinnedDataV[i] = (int) ceil(BndRankV[i] / Nr * NBin);
    AllBinnedDataV[i] = (int) ceil(AllRankV[i] / Nr * NBin);
    ExtraBinnedDataV[i] = (int) ceil(ExtraRankV[i] / Nr * NBin);
  }
  TFltV LklBayesFactorV; LklBayesFactorV.Gen(NBin);
  TFltV LklBayesDataV; LklBayesDataV.Gen(Nr);
  TFltV DnsBayesFactorV; DnsBayesFactorV.Gen(NBin);
  TFltV DnsBayesDataV; DnsBayesDataV.Gen(Nr);
  TFltV BndBayesFactorV; BndBayesFactorV.Gen(NBin);
  TFltV BndBayesDataV; BndBayesDataV.Gen(Nr);
  TFltV AllBayesFactorV; AllBayesFactorV.Gen(NBin);
  TFltV AllBayesDataV; AllBayesDataV.Gen(Nr);
  TFltV ExtraBayesFactorV; ExtraBayesFactorV.Gen(NBin);
  TFltV ExtraBayesDataV; ExtraBayesDataV.Gen(Nr);
  TFltV BayesDataSumV; BayesDataSumV.Gen(Nr);
  // estimated ranks, smaller is better (= closer to the top)
  // initialize the aggregated ranked list
  TFltV CRankLastV; CRankLastV.Gen(Nr);
  TFlt MaxBayes = -1e15;
  CRankV.Gen(Nr);
  for (int i = 0; i < Nr; i++) {
    CRankV[i] = (LklRankV[i] + DnsRankV[i] + BndRankV[i] + AllRankV[i]) / 4.;
    // Extra prioritization metric
    // CRankV[i] = (LklRankV[i] + DnsRankV[i] + BndRankV[i] + AllRankV[i] + ExtraRankV[i]) / 5.;
  }
  double CPrev = 0;
  TIntV OoV;
  for (int iter = 0; iter < MaxIter; iter++) {
    if (CorrStop - CPrev < 1e-15) {
      printf("Converged\n");
      break;
    }
    // assign the top np of aggregated predictions to be the positive class and the rest of the predictions to the negative class
    for (int i = 0; i < Nr; i++) { CRankLastV[i] = CRankV[i];}
     TCRankUtil::ArgSort(OoV, CRankV);
    for (int i = 0; i < Nr; i++) {
      if (i < Nrp) { CRankV[OoV[i]] = 1; }
      else { CRankV[OoV[i]] = 0; }
    }
    // Supervision
    //for (int i = 0; i < Nr; i++) {
    //  if (SideV[i] == 1.0) { CRankV[i] = 1;}
    //  else {CRankV[i] = 0;}
    //}
    // compute Bayes factors cumulatively starting from the top bin
    TCRank::GetBayesFactors(LklBayesFactorV, CRankV, LklBinnedDataV, NBin, Prior);
    TCRank::GetBayesFactors(DnsBayesFactorV, CRankV, DnsBinnedDataV, NBin, Prior);
    TCRank::GetBayesFactors(BndBayesFactorV, CRankV, BndBinnedDataV, NBin, Prior);
    TCRank::GetBayesFactors(AllBayesFactorV, CRankV, AllBinnedDataV, NBin, Prior);
    TCRank::GetBayesFactors(ExtraBayesFactorV, CRankV, ExtraBinnedDataV, NBin, Prior);
    // smooth using Tukey's running median
    // enforce a monotone decrease of Bayes factors
    TCRank::SmoothTukey(LklBayesFactorV);
    TCRank::SmoothTukey(DnsBayesFactorV);
    TCRank::SmoothTukey(BndBayesFactorV);
    TCRank::SmoothTukey(AllBayesFactorV);
    TCRank::SmoothTukey(ExtraBayesFactorV);
    LklBayesFactorV.Reverse(); TCRankUtil::CumMax(LklBayesFactorV); LklBayesFactorV.Reverse();
    DnsBayesFactorV.Reverse(); TCRankUtil::CumMax(DnsBayesFactorV); DnsBayesFactorV.Reverse();
    BndBayesFactorV.Reverse(); TCRankUtil::CumMax(BndBayesFactorV); BndBayesFactorV.Reverse();
    AllBayesFactorV.Reverse(); TCRankUtil::CumMax(AllBayesFactorV); AllBayesFactorV.Reverse();
    ExtraBayesFactorV.Reverse(); TCRankUtil::CumMax(ExtraBayesFactorV); ExtraBayesFactorV.Reverse();
    // if j-th entry in i-th ranking falls into k-th bin,
    // then bayes_data[j, i] should be the bayes factor of k-th bin in i-th ranking
    for (int i = 0; i < Nr; i++) {
      LklBayesDataV[i] = LklBayesFactorV[LklBinnedDataV[i]-1];
      DnsBayesDataV[i] = DnsBayesFactorV[DnsBinnedDataV[i]-1];
      BndBayesDataV[i] = BndBayesFactorV[BndBinnedDataV[i]-1];
      AllBayesDataV[i] = AllBayesFactorV[AllBinnedDataV[i]-1];
      ExtraBayesDataV[i] = ExtraBayesFactorV[ExtraBinnedDataV[i]-1];
      BayesDataSumV[i] = LklBayesDataV[i] + DnsBayesDataV[i] + BndBayesDataV[i] + AllBayesDataV[i];
      // Extra prioritization metric
      // BayesDataSumV[i] = LklBayesDataV[i] + DnsBayesDataV[i] + BndBayesDataV[i] + AllBayesDataV[i] + ExtraBayesDataV[i];
      //MaxBayes = fmax(MaxBayes, BayesDataSumV[i]);
      if (MaxBayes < BayesDataSumV[i]) {
        MaxBayes = BayesDataSumV[i];
      }
    }
    TCRankUtil::RankData(CRankV, BayesDataSumV);
    CPrev = TCRankUtil::PearsonCorrelation(CRankV, CRankLastV);
    printf("Correlation with previous iteration: %f\n", CPrev);
  }
  // return CRank scores for the final prioritized list; Bayes factor log space
  for (int i = 0; i < Nr; i++) {
    CRankV[i] = BayesDataSumV[i] / MaxBayes;
  }
}

/// smooth InV vector using Tukey's running median
void TCRank::SmoothTukey(TFltV& InV) {
  int nChanges = 0;
  TFltV TmpV;
  TFltV MV;
  while (nChanges > 0) {
    TmpV.Clr();
    for (int i = 0; i < InV.Len(); i++) { TmpV.Add(InV[i]); }
    for (int i = 1; i < InV.Len() - 1; i++) {
      MV.Clr();
      MV.Add(TmpV[i-1]); MV.Add(TmpV[i]); MV.Add(TmpV[i+1]);
      InV[i] = TMom(MV).GetMedian();
    }
    MV.Clr();
    MV.Add(TmpV[0]); MV.Add(InV[1]); MV.Add(3 * InV[1] - 2 * InV[2]);
    InV[0] = TMom(MV).GetMedian();
    MV.Clr();
    MV.Add(TmpV.Last()); MV.Add(InV.LastLast()); MV.Add(3 * InV.LastLast() - 2 * InV[InV.Len() - 3]);
    InV[InV.Len() - 1] = TMom(MV).GetMedian();
    nChanges = 0;
    for (int i = 0; i < InV.Len(); i++) {
        nChanges += InV[i] != TmpV[i] ? 1 : 0;
    }
  }
}

/// compute Bayes factors cumulatively starting from the top
void TCRank::GetBayesFactors(TFltV& BayesFactorV, const TFltV& RankV, const TIntV& BinnedDataV, const int& NBin, const double& Prior) {
  for (int bin = 1; bin < NBin + 1; bin++) {
    double Tpr = 0; double Fpr = 0;
    for (int i = 0; i < BinnedDataV.Len(); i++) {
      Tpr += BinnedDataV[i].Val <= bin ? RankV[i].Val : 0;
      Fpr += BinnedDataV[i].Val <= bin ? (1 - RankV[i].Val) : 0;
    }
    BayesFactorV[bin - 1] = log((Tpr + 1.) / (Fpr + 1.) / (Prior / (1. - Prior)));
  }
}

/// connect members of a given community by Erdos-Renyi
void TCRank::RndConnectInsideCommunity(PUNGraph& Graph, const TIntV& CmtyV, const double& Prob, TRnd& Rnd){
  int CNodes = CmtyV.Len(), CEdges;
  if (CNodes < 20) {
    CEdges = (int) Rnd.GetBinomialDev(Prob, CNodes * (CNodes-1) / 2);
  } else {
    CEdges = (int) (Prob * CNodes * (CNodes - 1) / 2);
  }
  THashSet<TIntPr> NewEdgeSet(CEdges);
  for (int edge = 0; edge < CEdges; ) {
    int SrcNId = CmtyV[Rnd.GetUniDevInt(CNodes)];
    int DstNId = CmtyV[Rnd.GetUniDevInt(CNodes)];
    if (SrcNId > DstNId) { Swap(SrcNId,DstNId); }
    if (SrcNId != DstNId && ! NewEdgeSet.IsKey(TIntPr(SrcNId, DstNId))) { // is new edge
      NewEdgeSet.AddKey(TIntPr(SrcNId, DstNId));
      Graph->AddEdge(SrcNId, DstNId);
      edge++; 
    } 
  }
}

void TCRank::ProbaConnectInsideCommunity(PUNGraph& Graph, const TIntV& CmtyV, const TVec<TFltV>& ProbaVV, TRnd& Rnd){
  THashSet<TIntPr> NewEdgeSet;
  for (int i = 0; i < ProbaVV.Len(); i++) {
    for (int j = 0; j < ProbaVV[i].Len(); j++) {
      if (Rnd.GetUniDev() < ProbaVV[i][j]) {
        int SrcNId = CmtyV[i];
        int DstNId = CmtyV[j];
        if (SrcNId > DstNId) { Swap(SrcNId,DstNId); }
        if (SrcNId != DstNId && ! NewEdgeSet.IsKey(TIntPr(SrcNId, DstNId))) {
          NewEdgeSet.AddKey(TIntPr(SrcNId, DstNId));
          Graph->AddEdge(SrcNId, DstNId);
        }
      }
    }
  }
}

PUNGraph TCRank::GenProbaModel(TVec<TIntV>& CmtyVV, const THash<TIntPr,TFlt>& CProbaH, const THash<TIntPr,TFlt>& EdgeProbaH, const THash<TIntTr,TFlt>& CEdgeProbaH, const PUNGraph& OrigGraph, TRnd& Rnd){
  PUNGraph G = TUNGraph::New(100 * CmtyVV.Len(), -1);
  for (TUNGraph::TNodeI NI = OrigGraph->BegNI(); NI < OrigGraph->EndNI(); NI++) {
    G->AddNode(NI.GetId());
  }
  // community-specific edge probabilities
  for (int i = 0; i < CmtyVV.Len(); i++) {
    TVec<TFltV> ProbaVV;
    for (int u = 0; u < CmtyVV[i].Len(); u++) {
      TFltV ProbaV;
      for (int v = 0; v < CmtyVV[i].Len(); v++) {
        TIntTr K1 = TIntTr(CmtyVV[i][u], CmtyVV[i][v], i);
        TIntPr K2 = TIntPr(CmtyVV[i][u], i);
        TIntPr K3 = TIntPr(CmtyVV[i][v], i);
        if (CEdgeProbaH.IsKey(K1) && CProbaH.IsKey(K2) && CProbaH.IsKey(K3)) {
          double CProba = CEdgeProbaH.GetDat(K1) * CProbaH.GetDat(K2) * CProbaH.GetDat(K3);
          ProbaV.Add(CProba);
        }
        else {
          ProbaV.Add(0);
        }
      }
      ProbaVV.Add(ProbaV);
    }
    ProbaConnectInsideCommunity(G, CmtyVV[i], ProbaVV, Rnd);
  }
  // edge probabilities
  for (THash<TIntPr,TFlt>::TIter HI = EdgeProbaH.BegI(); HI < EdgeProbaH.EndI(); HI++){
    double Proba = HI.GetDat();
    if (Rnd.GetUniDev() < Proba) {
      int SrcNId = HI.GetKey().Val1;
      int DstNId = HI.GetKey().Val2;
      if (SrcNId > DstNId) { Swap(SrcNId, DstNId); }
      if (SrcNId != DstNId && !G->IsEdge(SrcNId, DstNId)) {
        G->AddEdge(SrcNId, DstNId);
      }
    }
  }
  G->Defrag();
  return G;
}

PUNGraph TCRank::GenProbaModelAlpha(TVec<TIntV>& CmtyVV, const THash<TIntPr,TFlt>& CProbaH, const THash<TIntPr,TFlt>& EdgeProbaH, const THash<TIntTr,TFlt>& CEdgeProbaH, const PUNGraph& OrigGraph, const double Alpha, TRnd& Rnd){
  PUNGraph G = TUNGraph::New(100 * CmtyVV.Len(), -1);
  int TargetEdges = OrigGraph->GetEdges();
  for (TUNGraph::TNodeI NI = OrigGraph->BegNI(); NI < OrigGraph->EndNI(); NI++) {
    G->AddNode(NI.GetId());
  }
  // community-specific edge probabilities
  for (int i = 0; i < CmtyVV.Len(); i++) {
    TVec<TFltV> ProbaVV;
    for (int u = 0; u < CmtyVV[i].Len(); u++) {
      TFltV ProbaV;
      for (int v = 0; v < CmtyVV[i].Len(); v++) {
        TIntTr K1 = TIntTr(CmtyVV[i][u], CmtyVV[i][v], i);
        TIntPr K2 = TIntPr(CmtyVV[i][u], i);
        TIntPr K3 = TIntPr(CmtyVV[i][v], i);
        if (CEdgeProbaH.IsKey(K1) && CProbaH.IsKey(K2) && CProbaH.IsKey(K3)) {
          double CProba = CEdgeProbaH.GetDat(K1) * CProbaH.GetDat(K2) * CProbaH.GetDat(K3);
          TUNGraph::TNodeI UNI = OrigGraph->GetNI(CmtyVV[i][u]);
          double euv = UNI.GetOutDeg();
          TUNGraph::TNodeI VNI = OrigGraph->GetNI(CmtyVV[i][v]);
          euv *= VNI.GetOutDeg() / (2 * TargetEdges);
          double CProbaPert = CProba * (1 - Alpha) + (1 - CProba) * (1 - pow((1 - euv / TargetEdges), Alpha * TargetEdges));
          ProbaV.Add(CProbaPert);
        }
        else {
          ProbaV.Add(0);
        }
      }
      ProbaVV.Add(ProbaV);
    }
    ProbaConnectInsideCommunity(G, CmtyVV[i], ProbaVV, Rnd);
  }
  // edge probabilities
  for (THash<TIntPr,TFlt>::TIter HI = EdgeProbaH.BegI(); HI < EdgeProbaH.EndI(); HI++){
    int SrcNId = HI.GetKey().Val1;
    int DstNId = HI.GetKey().Val2;
    double Proba = HI.GetDat();
    TUNGraph::TNodeI SNI = OrigGraph->GetNI(SrcNId);
    double euv = SNI.GetOutDeg();
    TUNGraph::TNodeI DNI = OrigGraph->GetNI(DstNId);
    euv *= DNI.GetOutDeg() / (2 * TargetEdges);
    double ProbaPert = Proba * (1 - Alpha) + (1 - Proba) * (1 - pow((1 - euv / TargetEdges), Alpha * TargetEdges));
    if (Rnd.GetUniDev() < ProbaPert) {
      if (SrcNId > DstNId) { Swap(SrcNId, DstNId); }
      if (SrcNId != DstNId && !G->IsEdge(SrcNId, DstNId)) {
        G->AddEdge(SrcNId, DstNId);
      }
    }
  }
  // add remaining edges with a given degree sequence
  int NAdd = OrigGraph->GetEdges() - G->GetEdges();
  // printf("NAdd: %d\n", NAdd);
  if (NAdd > 0 && Alpha > 0) {
    TIntV DegSegV;
    for (int u = 0; u < OrigGraph->GetNodes(); u++) {
      int Deg = OrigGraph->GetNI(u).GetDeg();
      DegSegV.Add(Deg);
    }
    AddConfModel(G, DegSegV, NAdd, Rnd);
  }
  G->Defrag();
  return G;
}

void TCRank::GenCmtyProbV(TFltV& CProbV, TVec<TIntV>& CmtyVV, const double& DensityCoef, const double& ScaleCoef) {
  double Prob;
  for (int i = 0; i < CmtyVV.Len(); i++) {
    Prob = ScaleCoef * pow(double(CmtyVV[i].Len()), - DensityCoef);
    if (Prob > 1.0) { Prob = 1; }
    CProbV.Add(Prob);
  }
}

void TCRank::AddConfModel(PUNGraph& Graph, const TIntV& DegSeqV, const int& NAdd, TRnd& Rnd) {
  //printf("Conf model begins\n");
  const int Nodes = DegSeqV.Len();
  TIntV NIdDegV(DegSeqV.Len(), 0);
  for (int node = 0; node < Nodes; node++) {
    for (int d = 0; d < DegSeqV[node]; d++) {
      NIdDegV.Add(node);
    }
  }
  int edges=0;
  NIdDegV.Shuffle(Rnd);
  TIntPrSet EdgeH(NAdd); // set of all edges to be added
  int u=0, v=0;
  for (int c = 0; c < NAdd; c++) {
    u = Rnd.GetUniDevInt(NIdDegV.Len());
    while ((v = Rnd.GetUniDevInt(NIdDegV.Len())) == u) { }
    if (u > v) { Swap(u, v); }
    const int E1 = NIdDegV[u];
    const int E2 = NIdDegV[v];
    if (v == NIdDegV.Len()-1) { NIdDegV.DelLast(); }
    else { NIdDegV[v] = NIdDegV.Last();  NIdDegV.DelLast(); }
    if (u == NIdDegV.Len()-1) { NIdDegV.DelLast(); }
    else { NIdDegV[u] = NIdDegV.Last();  NIdDegV.DelLast(); }
    if (E1 == E2 || EdgeH.IsKey(TIntPr(E1, E2))) { continue; }
    EdgeH.AddKey(TIntPr(E1, E2));
    Graph->AddEdge(E1, E2);
    edges++;
    //if (c % 10000 == 0) { printf("Configuration model: iter %d: edges: %d, left: %d\n", c, edges,(NAdd-edges)); }
  }
  //printf("Conf model completed (%d nodes %d edges)\n", Graph->GetNodes(), Graph->GetEdges());
}

PUNGraph TCRank::GenAGM(PUNGraph& Graph, TVec<TIntV>& CmtyVV, TRnd& Rnd){
  // DensityCoef a: Power-law coefficient a of density (density ~ N^(-a)
  // ScaleCoef s: Scaling coefficient s of density (density ~ s)
  const double DensityCoef = 0.6;
  PUNGraph TryG = TCRank::GenAGM(CmtyVV, DensityCoef, 1.0, Rnd);
  double ScaleCoef = (double) Graph->GetEdges() / (double) TryG->GetEdges();
  TFltV CProbV;
  TCRank::GenCmtyProbV(CProbV, CmtyVV, DensityCoef, ScaleCoef);
  return TCRank::GenAGM(CmtyVV, CProbV, Rnd);
}

PUNGraph TCRank::GenAGM(TVec<TIntV>& CmtyVV, const double& DensityCoef, const double& ScaleCoef, TRnd& Rnd){
  TFltV CProbV;
  TCRank::GenCmtyProbV(CProbV, CmtyVV, DensityCoef, ScaleCoef);
  return TCRank::GenAGM(CmtyVV, CProbV, Rnd);
}

/// generate graph using the AGM model. CProbV = vector of Pc
PUNGraph TCRank::GenAGM(TVec<TIntV>& CmtyVV, const TFltV& CProbV, TRnd& Rnd, const double PNoCom){
  PUNGraph G = TUNGraph::New(100 * CmtyVV.Len(), -1);
  //printf("AGM begins\n");
  for (int i = 0; i < CmtyVV.Len(); i++) {
    TIntV& CmtyV = CmtyVV[i];
    for (int u = 0; u < CmtyV.Len(); u++) {
      if (G->IsNode(CmtyV[u])) { continue; }
      G->AddNode(CmtyV[u]);
    }
    double Prob = CProbV[i];
    RndConnectInsideCommunity(G, CmtyV, Prob, Rnd);
  }
  if (PNoCom > 0.0) { //if we want to connect nodes that do not share any community
    TIntSet NIDS;
    for (int c = 0; c < CmtyVV.Len(); c++) {
      for (int u = 0; u < CmtyVV[c].Len(); u++) {
        NIDS.AddKey(CmtyVV[c][u]);
      }
    }
    TIntV NIDV;
    NIDS.GetKeyV(NIDV);
    RndConnectInsideCommunity(G, NIDV, PNoCom, Rnd);
  }
  G->Defrag();
  //printf("AGM completed (%d nodes %d edges)\n", G->GetNodes(), G->GetEdges());
  return G;
}

PUNGraph TCRank::GenAGMAlpha(PUNGraph& Graph, TVec<TIntV >& CmtyVV, const double Alpha, TRnd& Rnd) {
  // DensityCoef a: Power-law coefficient a of density (density ~ N^(-a)
  // ScaleCoef s: Scaling coefficient s of density (density ~ s)
  const double DensityCoef = 0.6;
  PUNGraph TryG = TCRank::GenAGM(CmtyVV, DensityCoef, 1.0, Rnd);
  const double ScaleCoef = (double) Graph->GetEdges() / (double) TryG->GetEdges();
  TFltV CProbV;
  TCRank::GenCmtyProbV(CProbV, CmtyVV, DensityCoef, ScaleCoef);
  int Edges = Graph->GetEdges();
  TIntV DegSegV;
  for (int u = 0; u < Graph->GetNodes(); u++) {
    int Deg = Graph->GetNI(u).GetDeg();
    DegSegV.Add(Deg);
  }
  return TCRank::GenAGMAlpha(CmtyVV, CProbV, DegSegV, Edges, Alpha, Rnd);
}

/// generate graph using the AGM model. CProbV = vector of Pc, alpha = perturbation rate
PUNGraph TCRank::GenAGMAlpha(TVec<TIntV>& CmtyVV, const TFltV& CProbV, const TIntV& DegSegV, const int& Edges, const double& Alpha, TRnd& Rnd, const double PNoCom){
  PUNGraph G = TUNGraph::New(100 * CmtyVV.Len(), -1);
  //printf("Alpha: %f\n", Alpha);
  //printf("AGM begins\n");
  for (int i = 0; i < CmtyVV.Len(); i++) {
    TIntV& CmtyV = CmtyVV[i];
    for (int u = 0; u < CmtyV.Len(); u++) {
      if (G->IsNode(CmtyV[u])) { continue; }
      G->AddNode(CmtyV[u]);
    }
    // Alpha = perturbation rate
    double Prob = (1 - Alpha) * CProbV[i];
    RndConnectInsideCommunity(G, CmtyV, Prob, Rnd);
  }
  if (PNoCom > 0.0) { //if we want to connect nodes that do not share any community
    TIntSet NIDS;
    for (int c = 0; c < CmtyVV.Len(); c++) {
      for (int u = 0; u < CmtyVV[c].Len(); u++) {
        NIDS.AddKey(CmtyVV[c][u]);
      }
    }
    TIntV NIDV;
    NIDS.GetKeyV(NIDV);
    RndConnectInsideCommunity(G, NIDV, PNoCom, Rnd);
  }
  //printf("AGM completed (%d nodes %d edges)\n", G->GetNodes(), G->GetEdges());
  // add any new nodes to G (node ids are ints)
  for (int u = 0; u < DegSegV.Len(); u++){
      if (G->IsNode(u)) { continue; }
      G->AddNode(u);
  }
  // add NAdd random edges with a given degree sequence
  int NAdd = Edges - G->GetEdges();
  //printf("NAdd: %d\n", NAdd);
  if (NAdd > 0 && Alpha > 0) {
    AddConfModel(G, DegSegV, NAdd, Rnd);
  }
  G->Defrag();
  return G;
}

/////////////////////////////////////////////////
// CRank prioritization model

void TCRankMetric::GetLikelihoodCmty(TFltV& LklV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV) {
  LklV.Gen(CmtyVV.Len());
  double Lkl;
  for (int c = 0; c < CmtyVV.Len(); c++) {
    TIntSet NIdSet(CmtyVV[c].Len());
    for (int i = 0; i < CmtyVV[c].Len(); i++) { // edges inside
        NIdSet.AddKey(CmtyVV[c][i]);
    }
    Lkl = 1;
    for (int i = 0; i < CmtyVV[c].Len(); i++) {
      TUNGraph::TNodeI NI = Graph->GetNI(CmtyVV[c][i]);
      int NbrIn = 0;
      for (int j = 0; j < NI.GetOutDeg(); j++) {
        if (NIdSet.IsKey(NI.GetOutNId(j))) { NbrIn += 1; }
      }
      // Additive Laplace smoothing, alpha=1
      Lkl *= (NbrIn + 1) / (double) (NI.GetOutDeg() + 1);
    }
    LklV[c] = (CmtyVV[c].Len() - 1) * log(Lkl) - log(CmtyVV[c].Len() * (CmtyVV[c].Len() - 1));
  }
}

void TCRankMetric::GetDensityCmty(TFltV& DnsV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV) {
  DnsV.Gen(CmtyVV.Len());
  double Dns;
  for (int c = 0; c < CmtyVV.Len(); c++) {
    PUNGraph SG = TSnap::GetSubGraph(Graph, CmtyVV[c]);
    double Nodes = CmtyVV[c].Len();
    Dns = SG->GetEdges() / (Nodes * (Nodes - 1));
    DnsV[c] = Dns;
  }
}

void TCRankMetric::GetBoundaryCmty(TFltV& BndV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV) {
  BndV.Gen(CmtyVV.Len());
  double Bnd;
  for (int c = 0; c < CmtyVV.Len(); c++) {
    TIntSet NIdSet(CmtyVV[c].Len());
    for (int i = 0; i < CmtyVV[c].Len(); i++) { // edges inside
        NIdSet.AddKey(CmtyVV[c][i]);
    }
    Bnd = 0;
    for (int i = 0; i < CmtyVV[c].Len(); i++) {
      TUNGraph::TNodeI NI = Graph->GetNI(CmtyVV[c][i]);
      int NbrOut = 0;
      for (int j = 0; j < NI.GetOutDeg(); j++) {
        if (! NIdSet.IsKey(NI.GetOutNId(j))) { NbrOut += 1; }
      }
      // Higher value is better
      Bnd -= NbrOut;
    }
    BndV[c] = Bnd / (double) CmtyVV[c].Len();
  }
}

void TCRankMetric::GetAllegianceCmty(TFltV& AllV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV) {
  AllV.Gen(CmtyVV.Len());
  double All;
  for (int c = 0; c < CmtyVV.Len(); c++) {
    TIntSet NIdSet(CmtyVV[c].Len());
    for (int i = 0; i < CmtyVV[c].Len(); i++) { // edges inside
        NIdSet.AddKey(CmtyVV[c][i]);
    }
    All = 0;
    for (int i = 0; i < CmtyVV[c].Len(); i++) {
      TUNGraph::TNodeI NI = Graph->GetNI(CmtyVV[c][i]);
      int NbrIn = 0, NbrOut = 0;
      for (int j = 0; j < NI.GetOutDeg(); j++) {
        if (NIdSet.IsKey(NI.GetOutNId(j))) { NbrIn += 1; }
        else { NbrOut += 1; }
      }
      All += NbrIn > NbrOut ? 1 : 0;
    }
    AllV[c] = All / (double) CmtyVV[c].Len();
  }
}

void TCRankMetric::GetPriorMetr(TFltV& MetricScoresV, const TFltV& OrigV, const TFltV& PertV) {
  MetricScoresV.Gen(OrigV.Len());
  double Tmp;
  for (int i = 0; i < OrigV.Len(); i++) {
    Tmp = OrigV[i] / (1. + fabs(OrigV[i] - PertV[i]));
    MetricScoresV[i] = Tmp;
  }
}

////////////////////////////////////////////////////////////////////////////////////
/// CRankUtil:: Utilities for CRank

/// compute cumulative max vector
void TCRankUtil::CumMax(TFltV& DatV) {
  double Mx = DatV[0];
  for (int i = 1; i < DatV.Len(); i++) {
    //DatV[i] = fmax(Mx, DatV[i].Val);
    if (DatV[i].Val < Mx) {
      DatV[i] = Mx;
    }
     Mx = DatV[i];
  }
}

double TCRankUtil::PearsonCorrelation(const TFltV& InA, const TFltV& InB) {
  TMom AM = TMom(InA);
  TMom BM = TMom(InB);
  double meanA = AM.GetMean(); double meanB = BM.GetMean();
  double stDevA = AM.GetSDev(); double stDevB = BM.GetSDev();
  double pearson = 0;
  for (int i = 0; i < InA.Len(); i++) {
    pearson += (InA[i] - meanA) * (InB[i] - meanB);
  }
  pearson /= InA.Len() * stDevA * stDevB;
  return pearson;
}

void TCRankUtil::ArgSort(TIntV& OutV, const TFltV InV) {
  OutV.Gen(InV.Len());
  TFltV TmpV;
  for (int i = 0; i < InV.Len(); i++) {
    TmpV.Add(InV[i]);
  }
  TmpV.Sort();
  THash<TInt,TFlt> MapOutH;
  for (int i = 0; i < TmpV.Len(); i++) { MapOutH.AddDat(i, TmpV[i]); }
  THash<TFlt,TIntV> MapInVH;
  for (int i = 0; i < InV.Len(); i++) {
    if (!MapInVH.IsKey(InV[i])) {
      MapInVH.AddDat(InV[i], TIntV());
    }
    MapInVH.GetDat(InV[i]).Add(i);
  }
  for (TInt i = 0; i < InV.Len(); i++) {
    TFlt k = MapOutH.GetDat(i);
    OutV[i] = MapInVH.GetDat(k)[0];
    MapInVH.GetDat(k).Del(0);
    // printf("%d\n", OutV[i].Val);
  }
}

/// assign ranks to data, dealing with ties appropriately, largest value gets the smallest rank, ranks begin at 1 / n.
void TCRankUtil::RankData(TFltV& OutV, const TFltV InV) {
  for (int i = 0; i < InV.Len(); i++) {
    OutV[i] = InV[i];
  }
  OutV.Sort(false);
  THash<TFlt,TFlt> MapH;
  THash<TFlt,TInt> CountH;
  for (int i = 0; i < OutV.Len(); i++) {
    MapH.AddDat(OutV[i], i + 1);
    if (CountH.IsKey(OutV[i])) { CountH.AddDat(OutV[i], CountH.GetDat(OutV[i]) + 1); }
    else { CountH.AddDat(OutV[i], 1); }
  }
  // the average of the ranks that would have been assigned to all the tied values is assigned to each value
  for (THash<TFlt,TFlt>::TIter HI = MapH.BegI(); HI < MapH.EndI(); HI++) {
    TFlt Oval = HI.GetKey();
    if (CountH.GetDat(Oval) > 1) {
      // printf("%f\n", MapH.GetDat(Oval).Val);
      // printf("%d\n", CountH.GetDat(Oval).Val);
      TFlt AvgRank = (2 * MapH.GetDat(Oval) - CountH.GetDat(Oval) + 1) / 2.;
      MapH.AddDat(Oval, AvgRank);
    }
  }
  for (int i = 0; i < InV.Len(); i++) {
    OutV[i] = MapH.GetDat(InV[i]);
    // printf("%f -> %f\n", InV[i].Val, OutV[i].Val);
  }
}

/// load bipartite community affiliation graph from text file (each row contains the member node IDs for each community)
void TCRankUtil::LoadCmtyVV(const TStr& InFNm, TVec<TIntV>& CmtyVV) {
  CmtyVV.Gen(Kilo(100), 0);
  TSsParser Ss(InFNm, ssfWhiteSep);
  while (Ss.Next()) {
    if(Ss.GetFlds() > 0) {
      TIntV CmtyV;
      for (int i = 0; i < Ss.GetFlds(); i++) {
        if (Ss.IsInt(i)) {
          CmtyV.Add(Ss.GetInt(i));
        }
      }
      CmtyVV.Add(CmtyV);
    }
  }
  CmtyVV.Pack();
  printf("community loading completed (%d communities)\n",CmtyVV.Len());

}

/// load bipartite community affiliation graph from text file (each row contains the member node IDs for each community)
void TCRankUtil::LoadCmtyVV(const TStr& InFNm, TVec<TIntV>& CmtyVV, TStrHash<TInt>& StrToNIdH, const int BeginCol, const int MinSz, const TSsFmt Sep) {
  CmtyVV.Gen(Kilo(100), 0);
  TSsParser Ss(InFNm, Sep);
  while (Ss.Next()) {
    if(Ss.GetFlds() > BeginCol) {
      TIntV CmtyV;
      for (int i = BeginCol; i < Ss.GetFlds(); i++) {
        if (StrToNIdH.IsKey(Ss.GetFld(i))) {
          CmtyV.Add(StrToNIdH.GetKeyId(Ss.GetFld(i)));
        }
      }
      if (CmtyV.Len() < MinSz) { continue; }
      CmtyVV.Add(CmtyV);
    }
  }
  CmtyVV.Pack();
  printf("community loading completed (%d communities)\n",CmtyVV.Len());
}

/// dump bipartite community affiliation into a text file
void TCRankUtil::DumpCmtyVV(const TStr& OutFNm, const TVec<TIntV>& CmtyVV) {
  FILE* F = fopen(OutFNm.CStr(),"wt");
  for (int i = 0; i < CmtyVV.Len(); i++) {
    for (int j = 0; j < CmtyVV[i].Len(); j++) {
      fprintf(F,"%d\t", (int) CmtyVV[i][j]);
    }
    fprintf(F,"\n");
  }
  fclose(F);
}

/// dump bipartite community affiliation into a text file with node names
void TCRankUtil::DumpCmtyVV(const TStr OutFNm, TVec<TIntV>& CmtyVV, TIntStrH& NIDNmH) {
  FILE* F = fopen(OutFNm.CStr(), "wt");
  for (int c = 0; c < CmtyVV.Len(); c++) {
    for (int u = 0; u < CmtyVV[c].Len(); u++) {
      if (NIDNmH.IsKey(CmtyVV[c][u])){
        fprintf(F, "%s\t", NIDNmH.GetDat(CmtyVV[c][u]).CStr());
      }
      else {
        fprintf(F, "%d\t", (int) CmtyVV[c][u]);
      }
    }
    fprintf(F, "\n");
  }
  fclose(F);
}

/// total number of memberships (== sum of the sizes of communities)
int TCRankUtil::TotalMemberships(const TVec<TIntV>& CmtyVV){
  int M = 0;
  for (int i = 0; i < CmtyVV.Len(); i++) {
    M += CmtyVV[i].Len();
  }
  return M;
}

/// get hash table of <Node ID, membership size>
void TCRankUtil::GetNodeMembership(TIntH& NIDComVH, const THash<TInt,TIntV >& CmtyVH) {
  NIDComVH.Clr();
  for (THash<TInt,TIntV>::TIter HI = CmtyVH.BegI(); HI < CmtyVH.EndI(); HI++){
    for (int j = 0;j < HI.GetDat().Len(); j++) {
      int NID = HI.GetDat()[j];
      NIDComVH.AddDat(NID)++;
    }
  }
}

/// get hash table of <Node ID, community IDs which node belongs to>
void TCRankUtil::GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const TVec<TIntV>& CmtyVV) {
  NIDComVH.Clr();
  for (int i = 0; i < CmtyVV.Len(); i++){
    int CID = i;
    for (int j = 0; j < CmtyVV[i].Len(); j++) {
      int NID = CmtyVV[i][j];
      NIDComVH.AddDat(NID).AddKey(CID);
    }
  }
}

/// get hash table of <Node ID, community IDs which node belongs to>. Some nodes in NIDV might belong to no community
void TCRankUtil::GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const TVec<TIntV>& CmtyVV, const TIntV& NIDV) {
  NIDComVH.Clr();
  for (int u = 0; u < NIDV.Len(); u++) {
    NIDComVH.AddDat(NIDV[u]);
  }
  for (int i = 0; i < CmtyVV.Len(); i++){
    int CID = i;
    for (int j = 0; j < CmtyVV[i].Len(); j++) {
      int NID = CmtyVV[i][j];
      NIDComVH.AddDat(NID).AddKey(CID);
    }
  }
}

void TCRankUtil::GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const TVec<TIntSet>& CmtyVV) {
  for (int i = 0; i < CmtyVV.Len(); i++){
    int CID = i;
    for (TIntSet::TIter SI = CmtyVV[i].BegI(); SI < CmtyVV[i].EndI(); SI++) {
      int NID = SI.GetKey();
      NIDComVH.AddDat(NID).AddKey(CID);
    }
  }
}

void TCRankUtil::GetNodeMembership(THash<TInt,TIntSet >& NIDComVH, const THash<TInt,TIntV>& CmtyVH) {
  for (THash<TInt,TIntV>::TIter HI = CmtyVH.BegI(); HI < CmtyVH.EndI(); HI++){
    int CID = HI.GetKey();
    for (int j = 0; j < HI.GetDat().Len(); j++) {
      int NID = HI.GetDat()[j];
      NIDComVH.AddDat(NID).AddKey(CID);
    }
  }
}

void TCRankUtil::GetNodeMembership(THash<TInt,TIntV>& NIDComVH, const THash<TInt,TIntV>& CmtyVH) {
  for (int i = 0; i < CmtyVH.Len(); i++){
    int CID = CmtyVH.GetKey(i);
    for (int j = 0; j < CmtyVH[i].Len(); j++) {
      int NID = CmtyVH[i][j];
      NIDComVH.AddDat(NID).Add(CID);
    }
  }
}

void TCRankUtil::GetNodeMembership(THash<TInt,TIntV >& NIDComVH, const TVec<TIntV>& CmtyVV) {
  THash<TInt,TIntV> CmtyVH;
  for (int i = 0; i < CmtyVV.Len(); i++) {
    CmtyVH.AddDat(i, CmtyVV[i]);
  }
  GetNodeMembership(NIDComVH, CmtyVH);
}

int TCRankUtil::Intersection(const TIntV& C1, const TIntV& C2) {
  TIntSet S1(C1), S2(C2);
  return TCRankUtil::Intersection(S1, S2);
}

void TCRankUtil::GetIntersection(const THashSet<TInt>& A, const THashSet<TInt>& B, THashSet<TInt>& C) {
  C.Gen(A.Len());
  if (A.Len() < B.Len()) {
    for (THashSetKeyI<TInt> it = A.BegI(); it < A.EndI(); it++) 
      if (B.IsKey(it.GetKey())) C.AddKey(it.GetKey());
  } else {
    for (THashSetKeyI<TInt> it = B.BegI(); it < B.EndI(); it++) 
      if (A.IsKey(it.GetKey())) C.AddKey(it.GetKey());
  }
}

int TCRankUtil::Intersection(const THashSet<TInt>& A, const THashSet<TInt>& B) {
  int n = 0;
  if (A.Len() < B.Len()) {
    for (THashSetKeyI<TInt> it = A.BegI(); it < A.EndI(); it++) 
      if (B.IsKey(it.GetKey())) n++;
  } else {
    for (THashSetKeyI<TInt> it = B.BegI(); it < B.EndI(); it++) 
      if (A.IsKey(it.GetKey())) n++;
  }
  return n;
}

/// save graph into a gexf file which Gephi can read
void TCRankUtil::SaveGephi(const TStr& OutFNm, const PUNGraph& G, const TVec<TIntV>& CmtyVVAtr, const double MaxSz, const double MinSz, const TIntStrH& NIDNameH, const THash<TInt, TIntTr>& NIDColorH ) {
  THash<TInt,TIntV> NIDComVHAtr;
  TCRankUtil::GetNodeMembership(NIDComVHAtr, CmtyVVAtr);

  FILE* F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(F, "<gexf xmlns='http://www.gexf.net/1.2draft' xmlns:viz='http://www.gexf.net/1.1draft/viz' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xsi:schemaLocation='http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd' version='1.2'>\n");
  fprintf(F, "\t<graph mode='static' defaultedgetype='undirected'>\n");
  if (CmtyVVAtr.Len() > 0) {
    fprintf(F, "\t<attributes class='node'>\n");
    for (int c = 0; c < CmtyVVAtr.Len(); c++) {
      fprintf(F, "\t\t<attribute id='%d' title='c%d' type='boolean'>", c, c);
      fprintf(F, "\t\t<default>false</default>\n");
      fprintf(F, "\t\t</attribute>\n");
    }
    fprintf(F, "\t</attributes>\n");
  }
  fprintf(F, "\t\t<nodes>\n");
  for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    int NID = NI.GetId();
    TStr Label = NIDNameH.IsKey(NID)? NIDNameH.GetDat(NID): "";
    TIntTr Color = NIDColorH.IsKey(NID)? NIDColorH.GetDat(NID) : TIntTr(120, 120, 120);

    double Size = MinSz;
    double SizeStep = (MaxSz - MinSz) / (double) CmtyVVAtr.Len();
    if (NIDComVHAtr.IsKey(NID)) {
      Size = MinSz +  SizeStep *  (double) NIDComVHAtr.GetDat(NID).Len();
    }
    double Alpha = 1.0;
    fprintf(F, "\t\t\t<node id='%d' label='%s'>\n", NID, Label.CStr());
    fprintf(F, "\t\t\t\t<viz:color r='%d' g='%d' b='%d' a='%.1f'/>\n", Color.Val1.Val, Color.Val2.Val, Color.Val3.Val, Alpha);
    fprintf(F, "\t\t\t\t<viz:size value='%.3f'/>\n", Size);
    //specify attributes
    if (NIDComVHAtr.IsKey(NID)) {
      fprintf(F, "\t\t\t\t<attvalues>\n");
      for (int c = 0; c < NIDComVHAtr.GetDat(NID).Len(); c++) {
        int CID = NIDComVHAtr.GetDat(NID)[c];
        fprintf(F, "\t\t\t\t\t<attvalue for='%d' value='true'/>\n", CID);
      }
      fprintf(F, "\t\t\t\t</attvalues>\n");
    }

    fprintf(F, "\t\t\t</node>\n");
  }
  fprintf(F, "\t\t</nodes>\n");
  //plot edges
  int EID = 0;
  fprintf(F, "\t\t<edges>\n");
  for (TUNGraph::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
    fprintf(F, "\t\t\t<edge id='%d' source='%d' target='%d'/>\n", EID++, EI.GetSrcNId(), EI.GetDstNId());
  }
  fprintf(F, "\t\t</edges>\n");
  fprintf(F, "\t</graph>\n");
  fprintf(F, "</gexf>\n");
  fclose(F);
}

/// Save bipartite community affiliation into gexf file
void TCRankUtil::SaveBipartiteGephi(const TStr& OutFNm, const TIntV& NIDV, const TVec<TIntV>& CmtyVV, const double MaxSz, const double MinSz, const TIntStrH& NIDNameH, const THash<TInt, TIntTr>& NIDColorH, const THash<TInt, TIntTr>& CIDColorH ) {
  /// Plot bipartite graph
  if (CmtyVV.Len() == 0) { return; }
  double NXMin = 0.1, YMin = 0.1, NXMax = 250.00, YMax = 30.0;
  double CXMin = 0.3 * NXMax, CXMax = 0.7 * NXMax;
  double CStep = (CXMax - CXMin) / (double) CmtyVV.Len(), NStep = (NXMax - NXMin) / (double) NIDV.Len();
  THash<TInt,TIntV> NIDComVH;
  TCRankUtil::GetNodeMembership(NIDComVH, CmtyVV);

  FILE* F = fopen(OutFNm.CStr(), "wt");
  fprintf(F, "<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(F, "<gexf xmlns='http://www.gexf.net/1.2draft' xmlns:viz='http://www.gexf.net/1.1draft/viz' xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance' xsi:schemaLocation='http://www.gexf.net/1.2draft http://www.gexf.net/1.2draft/gexf.xsd' version='1.2'>\n");
  fprintf(F, "\t<graph mode='static' defaultedgetype='directed'>\n");
  fprintf(F, "\t\t<nodes>\n");
  for (int c = 0; c < CmtyVV.Len(); c++) {
    int CID = c;
    double XPos = c * CStep + CXMin;
    TIntTr Color = CIDColorH.IsKey(CID)? CIDColorH.GetDat(CID) : TIntTr(120, 120, 120);
    fprintf(F, "\t\t\t<node id='C%d' label='C%d'>\n", CID, CID);
    fprintf(F, "\t\t\t\t<viz:color r='%d' g='%d' b='%d'/>\n", Color.Val1.Val, Color.Val2.Val, Color.Val3.Val);
    fprintf(F, "\t\t\t\t<viz:size value='%.3f'/>\n", MaxSz);
    fprintf(F, "\t\t\t\t<viz:shape value='square'/>\n");
    fprintf(F, "\t\t\t\t<viz:position x='%f' y='%f' z='0.0'/>\n", XPos, YMax); 
    fprintf(F, "\t\t\t</node>\n");
  }

  for (int u = 0;u < NIDV.Len(); u++) {
    int NID = NIDV[u];
    TStr Label = NIDNameH.IsKey(NID)? NIDNameH.GetDat(NID): "";
    double Size = MinSz;
    double XPos = NXMin + u * NStep;
    TIntTr Color = NIDColorH.IsKey(NID)? NIDColorH.GetDat(NID) : TIntTr(120, 120, 120);
    double Alpha = 1.0;
    fprintf(F, "\t\t\t<node id='%d' label='%s'>\n", NID, Label.CStr());
    fprintf(F, "\t\t\t\t<viz:color r='%d' g='%d' b='%d' a='%.1f'/>\n", Color.Val1.Val, Color.Val2.Val, Color.Val3.Val, Alpha);
    fprintf(F, "\t\t\t\t<viz:size value='%.3f'/>\n", Size);
    fprintf(F, "\t\t\t\t<viz:shape value='square'/>\n");
    fprintf(F, "\t\t\t\t<viz:position x='%f' y='%f' z='0.0'/>\n", XPos, YMin); 
    fprintf(F, "\t\t\t</node>\n");
  }
  fprintf(F, "\t\t</nodes>\n");
  fprintf(F, "\t\t<edges>\n");
  int EID = 0;
  for (int u = 0;u < NIDV.Len(); u++) {
    int NID = NIDV[u];
    if (NIDComVH.IsKey(NID)) {
      for (int c = 0; c < NIDComVH.GetDat(NID).Len(); c++) {
        int CID = NIDComVH.GetDat(NID)[c];
        fprintf(F, "\t\t\t<edge id='%d' source='C%d' target='%d'/>\n", EID++, CID, NID);
      }
    }
  }
  fprintf(F, "\t\t</edges>\n");
  fprintf(F, "\t</graph>\n");
  fprintf(F, "</gexf>\n");
}

/// Generate conductance scores for a set of communities
void TCRankUtil::GetConductanceCmty(TFltV& PhiV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV) {
  PhiV.Gen(CmtyVV.Len());
  double Phi;
  for (int c = 0; c < CmtyVV.Len(); c++) {
    TIntSet CmtyS;
    for (int u = 0; u < CmtyVV[c].Len(); u++) {
      CmtyS.AddKey(CmtyVV[c][u]);
    }
    Phi = GetConductance(Graph, CmtyS, Graph->GetEdges());
    PhiV[c] = Phi;
  }
}

double TCRankUtil::GetConductance(const PUNGraph& Graph, const TIntSet& CmtyS, const int Edges) {
  const int Edges2 = Edges >= 0 ? 2*Edges : Graph->GetEdges();
  int Vol = 0,  Cut = 0; 
  double Phi = 0.0;
  for (int i = 0; i < CmtyS.Len(); i++) {
    if (! Graph->IsNode(CmtyS[i])) { continue; }
    TUNGraph::TNodeI NI = Graph->GetNI(CmtyS[i]);
    for (int e = 0; e < NI.GetOutDeg(); e++) {
      if (! CmtyS.IsKey(NI.GetOutNId(e))) { Cut += 1; }
    }
    Vol += NI.GetOutDeg();
  }
  // get conductance
  if (Vol != Edges2) {
    if (2 * Vol > Edges2) { Phi = Cut / double (Edges2 - Vol); }
    else if (Vol == 0) { Phi = 0.0; }
    else { Phi = Cut / double(Vol); }
  } else {
    if (Vol == Edges2) { Phi = 1.0; }
  }
  return Phi;
}

/// Generate modularity scores for a set of communities
void TCRankUtil::GetModularityCmty(TFltV& ModV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV) {
  ModV.Gen(CmtyVV.Len());
  double Mod;
  for (int c = 0; c < CmtyVV.Len(); c++) {
    Mod = GetModularity(Graph, CmtyVV[c], Graph->GetEdges());
    ModV[c] = Mod;
  }
}

double TCRankUtil::GetModularity(const PUNGraph& Graph, const TIntV& NIdV, const int Edges) {
  double EdgesIn = 0.0, EEdgesIn = 0.0; // EdgesIn=2*number of edges inside the cluster, EEdgesIn=expected edges inside
  TIntSet NIdSet(NIdV.Len());
  for (int e = 0; e < NIdV.Len(); e++) { // edges inside
    NIdSet.AddKey(NIdV[e]);
  }
  for (int e1 = 0; e1 < NIdV.Len(); e1++) {
    TUNGraph::TNodeI NI = Graph->GetNI(NIdV[e1]);
    EEdgesIn += NI.GetOutDeg();
    for (int i = 0; i < NI.GetOutDeg(); i++) {
      if (NIdSet.IsKey(NI.GetOutNId(i))) { EdgesIn += 1; }
    }
  }
  EEdgesIn = EEdgesIn*EEdgesIn / (2.0*Edges);
  if ((EdgesIn - EEdgesIn) == 0) { return 0; }
  else { return (EdgesIn - EEdgesIn) / (2.0*Edges); } // modularity
}

void TCRankUtil::GetNbhCom(const PUNGraph& Graph, const int NID, TIntSet& NBCmtyS) {
  TUNGraph::TNodeI NI = Graph->GetNI(NID);
  NBCmtyS.Gen(NI.GetDeg());
  NBCmtyS.AddKey(NID);
  for (int e = 0; e < NI.GetDeg(); e++) {
    NBCmtyS.AddKey(NI.GetNbrNId(e));
  }
}

/// Generate random scores for a set of communities
void TCRankUtil::GetRandomCmty(TFltV& RndV, const PUNGraph& Graph, const TVec<TIntV>& CmtyVV, TRnd& Rnd) {
  RndV.Gen(CmtyVV.Len(), 0);
  for (int c = 0; c < CmtyVV.Len(); c++) {
    RndV.Add(Rnd.GetUniDev());
  }
}
