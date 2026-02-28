interface SiteNavProps {
  currentPage: 'landing' | 'visit';
  onNavigateLanding: () => void;
  onNavigateVisit: () => void;
  explorerUrl?: string;
  onExploreData?: () => void;
}

export default function SiteNav({
  currentPage,
  onNavigateLanding,
  onNavigateVisit,
  explorerUrl,
  onExploreData,
}: SiteNavProps) {
  const handleExplore = () => {
    if (explorerUrl) {
      window.open(explorerUrl, '_blank', 'noopener');
    } else if (onExploreData) {
      onExploreData();
    }
  };

  return (
    <nav className="chrome-header sticky top-0 z-50 px-6 py-3 flex items-center justify-between">
      <button
        onClick={onNavigateLanding}
        className="text-lg font-bold tracking-tight"
        style={{ color: 'var(--accent-1)', background: 'none', border: 'none', cursor: 'pointer' }}
      >
        San José
      </button>

      <div className="flex items-center gap-4">
        <button
          onClick={onNavigateLanding}
          className={`site-nav-link ${currentPage === 'landing' ? 'site-nav-link--active' : ''}`}
        >
          Home
        </button>
        <button
          onClick={onNavigateVisit}
          className={`site-nav-link ${currentPage === 'visit' ? 'site-nav-link--active' : ''}`}
        >
          Project Summary
        </button>
        {(explorerUrl || onExploreData) && (
          <button onClick={handleExplore} className="btn-primary">
            Explore Data
          </button>
        )}
      </div>
    </nav>
  );
}
