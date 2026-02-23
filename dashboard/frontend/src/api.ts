import axios from 'axios';

// Use direct backend URL for static file serving, or /api for dev server
const API_BASE = import.meta.env.DEV ? '/api' : 'http://localhost:8000';

const api = axios.create({
  baseURL: API_BASE,
});

export interface PrecursorSummary {
  precursor_id: number;
  raw_file: string;
  mz: number;
  charge: number;
  rt_seconds: number;
  mobility: number;
  ccs: number | null;
  n_engines: number;
  confidence_weight: number | null;
  sage_cosine: number | null;
  // FragPipe (10 columns)
  fragpipe_peptide: string | null;
  fragpipe_modified: string | null;
  fragpipe_protein: string | null;
  fragpipe_probability: number | null;
  fragpipe_pep: number | null;
  fragpipe_hyperscore: number | null;
  fragpipe_qvalue: number | null;
  fragpipe_rt: number | null;
  fragpipe_mz: number | null;
  fragpipe_mobility: number | null;
  // Sage (12 columns)
  sage_peptide: string | null;
  sage_modified: string | null;
  sage_protein: string | null;
  sage_qvalue: number | null;
  sage_pep: number | null;
  sage_hyperscore: number | null;
  sage_peptide_qvalue: number | null;
  sage_protein_qvalue: number | null;
  sage_rt: number | null;
  sage_mz: number | null;
  sage_mobility: number | null;
  sage_match_tier: string | null;
  // DIA-NN (13 columns)
  diann_peptide: string | null;
  diann_modified: string | null;
  diann_protein: string | null;
  diann_qvalue: number | null;
  diann_pep: number | null;
  diann_global_qvalue: number | null;
  diann_pg_qvalue: number | null;
  diann_rt: number | null;
  diann_mz: number | null;
  diann_mobility: number | null;
  diann_match_tier: string | null;
  diann_match_score: number | null;
  // Raw
  collision_energy: number | null;
  raw_intensity_meta: number | null;
  isolation_mz: number | null;
  // Quality metrics
  ms1_rt_sigma: number | null;
  ms1_rt_r2: number | null;
  ms1_im_sigma: number | null;
  ms1_im_r2: number | null;
  isotope_cosim: number | null;
}

export interface RawFileInfo {
  name: string;
  count: number;
}

export interface SageMatchedFragment {
  fragment_type: string;  // "b" or "y"
  ion_number: number;     // fragment_ordinals (1, 2, 3...)
  charge: number;         // fragment_charge
  mz_experimental: number;
  mz_calculated: number;
  intensity: number;
}

export interface PrecursorDetail extends PrecursorSummary {
  sage_modified: string | null;
  fragment_mz: number[];
  fragment_intensity: number[];
  fragment_mobility: number[];
  fragment_scan: number[];
  sage_matched_fragments: SageMatchedFragment[] | null;
  sage_cosine: number | null;
  xic_rt: number[];
  xic_intensity: number[];
  mobilogram_im: number[];
  mobilogram_intensity: number[];
  isotope_mz: number[];
  isotope_intensity: number[];
  raw_rt: number[];
  raw_mz: number[];
  raw_mobility: number[];
  raw_intensity: number[];
}

export interface StoreInfo {
  path: string;
  num_precursors: number;
  num_row_groups: number;
  columns: string[];
}

export interface StoreStats {
  total_precursors: number;
  by_engines: Record<string, number>;
  by_charge: Record<string, number>;
  mz_range: [number, number];
  raw_files: Record<string, number>;
}

// Collection mode types
export interface DatasetQuality {
  rt_r2_median: number | null;
  im_r2_median: number | null;
  pct_high_quality: number | null;
}

export interface DatasetSummary {
  accession: string;
  version: string;
  study_id: string;
  path: string;
  n_precursors: number;
  n_all_three: number | null;
  n_at_least_two: number | null;
  quality: DatasetQuality | null;
  added_at: string | null;
}

export interface DatasetOverlapSummary {
  dataset_a: string;
  dataset_b: string;
  precursors_a: number;
  precursors_b: number;
  shared_precursors: number;
  unique_precursors_a: number;
  unique_precursors_b: number;
  precursor_jaccard: number;
  peptides_a: number;
  peptides_b: number;
  shared_peptides: number;
  unique_peptides_a: number;
  unique_peptides_b: number;
  peptide_jaccard: number;
}

export interface StudySummary {
  id: string;
  title: string;
  organism: string | null;
  publication: string | null;
  description: string | null;
  n_datasets: number;
  n_total_precursors: number;
  datasets: string[];
}

export interface CollectionInfo {
  version: string;
  updated_at: string;
  n_studies: number;
  n_datasets: number;
  n_total_precursors: number;
}

export interface EngineOverlapStats {
  n_fragpipe: number;
  n_diann: number;
  n_sage: number;
  n_all_three: number;
  n_fp_dn_only: number;
  n_fp_sg_only: number;
  n_dn_sg_only: number;
  n_fragpipe_only: number;
  n_diann_only: number;
  n_sage_only: number;
  n_union: number;
  n_at_least_two: number;
  three_way_rate: number;
  at_least_two_rate: number;
}

export interface QualitySummaryData {
  n_0_engines: number;
  pct_0_engines: number;
  n_1_engines: number;
  pct_1_engines: number;
  n_2_engines: number;
  pct_2_engines: number;
  n_3_engines: number;
  pct_3_engines: number;
  ms1_rt_r2_mean: number | null;
  ms1_rt_r2_median: number | null;
  ms1_im_r2_mean: number | null;
  ms1_im_r2_median: number | null;
  isotope_cosim_mean: number | null;
  isotope_cosim_median: number | null;
  n_high_quality: number;
  pct_high_quality: number;
}

export interface RawFileSummaryData {
  name: string;
  count: number;
}

export interface QualityDistBucket {
  lower: number;
  upper: number;
  count: number;
  pct: number;
}

export interface DatasetFullSummary {
  accession: string | null;
  pipeline_version: string | null;
  generated_at: string | null;
  n_total_precursors: number;
  n_per_engine: Record<string, number> | null;
  n_unidentified: number;
  quality_summary: QualitySummaryData | null;
  overlap_stats: EngineOverlapStats | null;
  by_charge: Record<string, number>;
  mz_range: [number, number] | null;
  rt_range_minutes: [number, number] | null;
  mobility_range: [number, number] | null;
  collision_energy_range: [number, number] | null;
  raw_files: RawFileSummaryData[];
  n_unique_peptides: Record<string, number> | null;
  quality_distributions: Record<string, QualityDistBucket[]> | null;
  organism: string | null;
  study_id: string | null;
  group_id: string | null;
}

export interface AppStatus {
  status: string;
  mode: 'single' | 'collection';
  store_loaded: boolean;
  active_dataset: string | null;
  collection_loaded: boolean;
}

export async function loadStore(path: string) {
  const response = await api.post('/load', null, { params: { path } });
  return response.data;
}

export async function getStoreInfo(): Promise<StoreInfo> {
  const response = await api.get('/info');
  return response.data;
}

export async function getStats(): Promise<StoreStats> {
  const response = await api.get('/stats');
  return response.data;
}

export async function listPrecursors(params: {
  offset?: number;
  limit?: number;
  min_engines?: number;
  max_engines?: number;
  charge?: number;
  raw_file?: string;
  has_ms1?: boolean;
  has_raw_data?: boolean;
  sort_by?: string;
  sort_desc?: boolean;
  overlap_mode?: 'shared' | 'unique_a';
  overlap_dataset?: string;
}): Promise<PrecursorSummary[]> {
  const response = await api.get('/precursors', { params });
  return response.data;
}

export async function getRawFiles(): Promise<RawFileInfo[]> {
  const response = await api.get('/raw_files');
  return response.data;
}

export async function getPrecursor(id: number, rawFile?: string): Promise<PrecursorDetail> {
  const params: Record<string, string> = {};
  if (rawFile) params.raw_file = rawFile;
  const response = await api.get(`/precursor/${id}`, { params });
  return response.data;
}

// Collection mode API functions
export async function getAppStatus(): Promise<AppStatus> {
  const response = await api.get('/');
  return response.data;
}

export async function getCollectionInfo(): Promise<CollectionInfo> {
  const response = await api.get('/collection');
  return response.data;
}

export async function getStudies(): Promise<StudySummary[]> {
  const response = await api.get('/studies');
  return response.data;
}

export async function getStudyDatasets(studyId: string): Promise<DatasetSummary[]> {
  const response = await api.get(`/studies/${studyId}/datasets`);
  return response.data;
}

export async function getDatasetInfo(accession: string): Promise<DatasetSummary> {
  const response = await api.get(`/datasets/${accession}/info`);
  return response.data;
}

export async function getDatasetOverlap(datasetA: string, datasetB: string): Promise<DatasetOverlapSummary> {
  const response = await api.get('/datasets/overlap', {
    params: {
      dataset_a: datasetA,
      dataset_b: datasetB,
    },
  });
  return response.data;
}

export async function loadDataset(accession: string): Promise<{ status: string; accession: string; path: string; num_precursors: number }> {
  const response = await api.post(`/datasets/${accession}/load`);
  return response.data;
}

export async function getDatasetFullSummary(): Promise<DatasetFullSummary> {
  const response = await api.get('/summary');
  return response.data;
}

export async function getActiveDataset(): Promise<{ accession: string | null; path: string | null; loaded: boolean }> {
  const response = await api.get('/datasets/active');
  return response.data;
}

export default api;
