import { useState } from 'react';
import ReactMarkdown from 'react-markdown';
import remarkGfm from 'remark-gfm';
import SiteNav from './SiteNav';

import pitchRaw from '@docs/SAN_JOSE_PITCH.md?raw';
import projectPlanRaw from '@root/CLAUDIUS-PROTEOMICS.md?raw';
import runnerArchRaw from '@docs/RUNNER_ARCHITECTURE.md?raw';
import rawDataArchRaw from '@docs/RAW_DATA_ARCHITECTURE.md?raw';
import outputSchemaRaw from '@docs/RUNNER_OUTPUT_SCHEMA.md?raw';
import preScaleRaw from '@docs/PRE_SCALE_DECISIONS.md?raw';
import remediationRaw from '@docs/REMEDIATION_CHECKLIST.md?raw';
import datasetDefRaw from '@docs/DATASET_DEFINITION.md?raw';

interface Section {
  id: string;
  title: string;
  subtitle: string;
  content: string;
}

function stripLeadingH1(md: string): string {
  return md.replace(/^\s*#\s+[^\n]+\n*/, '');
}

const sections: Section[] = [
  {
    id: 'vision',
    title: 'Vision & Pitch',
    subtitle: 'Why San Jos\u00e9 exists and what it aims to achieve',
    content: stripLeadingH1(pitchRaw),
  },
  {
    id: 'project-plan',
    title: 'Project Plan',
    subtitle: 'Full technical plan and milestones',
    content: stripLeadingH1(projectPlanRaw),
  },
  {
    id: 'runner-architecture',
    title: 'Runner Architecture',
    subtitle: '6-step processing pipeline design',
    content: stripLeadingH1(runnerArchRaw),
  },
  {
    id: 'raw-data',
    title: 'Raw Data Architecture',
    subtitle: '4D signal storage, Parquet schemas, and dashboard API',
    content: stripLeadingH1(rawDataArchRaw),
  },
  {
    id: 'output-schema',
    title: 'Output Schema',
    subtitle: 'Parquet schemas, manifest format, and runner outputs',
    content: stripLeadingH1(outputSchemaRaw),
  },
  {
    id: 'pre-scale',
    title: 'Pre-Scale Decisions',
    subtitle: 'Architectural choices made before scaling up',
    content: stripLeadingH1(preScaleRaw),
  },
  {
    id: 'remediation',
    title: 'Remediation Checklist',
    subtitle: 'Known issues and fixes to apply before production',
    content: stripLeadingH1(remediationRaw),
  },
  {
    id: 'dataset-definition',
    title: 'Dataset Definition',
    subtitle: 'Dataset semantics and edge cases',
    content: stripLeadingH1(datasetDefRaw),
  },
];

interface BlueprintPageProps {
  onBack: () => void;
  onNavigateVisit: () => void;
  explorerUrl?: string;
  onExploreData?: () => void;
}

export default function BlueprintPage({
  onBack,
  onNavigateVisit,
  explorerUrl,
  onExploreData,
}: BlueprintPageProps) {
  const [activeSectionId, setActiveSectionId] = useState(sections[0].id);
  const [tocOpen, setTocOpen] = useState(false);

  const activeIndex = sections.findIndex((s) => s.id === activeSectionId);
  const active = sections[activeIndex];
  const prev = activeIndex > 0 ? sections[activeIndex - 1] : null;
  const next = activeIndex < sections.length - 1 ? sections[activeIndex + 1] : null;

  const navigateTo = (id: string) => {
    setActiveSectionId(id);
    setTocOpen(false);
    window.scrollTo({ top: 0 });
  };

  return (
    <div className="app-shell min-h-screen flex flex-col blueprint-page">
      <SiteNav
        currentPage="blueprint"
        onNavigateLanding={onBack}
        onNavigateVisit={onNavigateVisit}
        onNavigateBlueprint={() => {}}
        explorerUrl={explorerUrl}
        onExploreData={onExploreData}
      />

      <div className="blueprint-layout flex-1 flex">
        {/* Mobile TOC toggle */}
        <button
          className="blueprint-toc-toggle"
          onClick={() => setTocOpen(!tocOpen)}
          aria-label="Toggle table of contents"
        >
          <span className="blueprint-toc-toggle-icon">{tocOpen ? '\u2715' : '\u2630'}</span>
        </button>

        {/* TOC sidebar */}
        <aside className={`blueprint-toc${tocOpen ? ' blueprint-toc--open' : ''}`}>
          <div className="blueprint-toc-inner">
            <p className="subtle-label mb-4">Blueprint</p>
            <nav className="blueprint-toc-nav">
              {sections.map((s, i) => (
                <button
                  key={s.id}
                  onClick={() => navigateTo(s.id)}
                  className={`blueprint-toc-item${activeSectionId === s.id ? ' blueprint-toc-item--active' : ''}`}
                >
                  <span className="blueprint-toc-index mono">{String(i + 1).padStart(2, '0')}</span>
                  <span className="blueprint-toc-label">{s.title}</span>
                </button>
              ))}
            </nav>
          </div>
        </aside>

        {/* Main content — single section at a time */}
        <main className="blueprint-content">
          <div className="blueprint-content-inner">
            <header className="blueprint-section-header reveal-up">
              <span className="blueprint-section-index mono">
                {String(activeIndex + 1).padStart(2, '0')}
              </span>
              <div>
                <h1 className="blueprint-section-title">{active.title}</h1>
                <p className="blueprint-section-subtitle">{active.subtitle}</p>
              </div>
            </header>

            <div className="blueprint-prose">
              <ReactMarkdown remarkPlugins={[remarkGfm]}>{active.content}</ReactMarkdown>
            </div>

            {/* Prev / Next navigation */}
            <nav className="blueprint-pager">
              {prev ? (
                <button className="blueprint-pager-btn" onClick={() => navigateTo(prev.id)}>
                  <span className="blueprint-pager-dir">&larr; Previous</span>
                  <span className="blueprint-pager-title">{prev.title}</span>
                </button>
              ) : (
                <span />
              )}
              {next ? (
                <button
                  className="blueprint-pager-btn blueprint-pager-btn--next"
                  onClick={() => navigateTo(next.id)}
                >
                  <span className="blueprint-pager-dir">Next &rarr;</span>
                  <span className="blueprint-pager-title">{next.title}</span>
                </button>
              ) : (
                <span />
              )}
            </nav>
          </div>
        </main>
      </div>
    </div>
  );
}
