/* Modern teaching cards: the observed statistic, its null distribution and the tail area. */
(function (global) {
  'use strict';
  const raw = String.raw;
  const handled = new Set(['z-statistic', 't-statistic', 'pearson']);
  const retained = new Map();
  let current = null;
  const clamp01 = x => Math.max(0, Math.min(1, x));
  const esc = x => String(x ?? '').replace(/[&<>"']/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
  const format = (x, places = 4) => Number.isFinite(x) ? x.toLocaleString('ru-RU', {maximumFractionDigits: places}) : '—';
  const probability = p => p < 0.000001 ? '< 0,000001' : format(p, 6);
  const source = (title, locator, url) => ({title, author: 'NIST/SEMATECH', year: 'электронный справочник', locator, url, type: 'modern'});
  const nistNormal = source('Normal Distribution', 'Раздел 1.3.6.6.1: стандартная нормальная модель и распределение среднего.', 'https://www.itl.nist.gov/div898/handbook/eda/section3/eda3661.htm');
  const nistMean = source('Confidence Limits for the Mean', 'Раздел 1.3.5.2: одновыборочный t-критерий, стандартная ошибка и доверительный интервал.', 'https://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm');
  const nistT = source("Critical Values of the Student's t Distribution", 'Раздел 1.3.6.7.2: степени свободы и односторонние/двусторонние пороги.', 'https://www.itl.nist.gov/div898/handbook/eda/section3/eda3672.htm');
  const nistChi = source('Chi-Square Goodness-of-Fit Test', 'Раздел 1.3.5.15: ожидаемые частоты, число степеней свободы и верхний хвост.', 'https://www.itl.nist.gov/div898/handbook/eda/section3/eda35f.htm');

  // Lanczos log-gamma; continued fractions/series for regularized beta and gamma.
  function logGamma(z) {
    const c = [676.5203681218851, -1259.1392167224028, 771.32342877765313, -176.61502916214059, 12.507343278686905, -.13857109526572012, 9.984369578019572e-6, 1.5056327351493116e-7];
    if (z < .5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z);
    z -= 1; let x = .99999999999980993;
    for (let i = 0; i < c.length; i++) x += c[i] / (z + i + 1);
    const t = z + c.length - .5;
    return .5 * Math.log(2 * Math.PI) + (z + .5) * Math.log(t) - t + Math.log(x);
  }
  function betaFraction(a, b, x) {
    const tiny = 1e-300;
    let c = 1, d = 1 - (a + b) * x / (a + 1);
    if (Math.abs(d) < tiny) d = tiny;
    d = 1 / d; let h = d;
    for (let m = 1; m <= 500; m++) {
      let aa = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m));
      d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny;
      c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1));
      d = 1 + aa * d; if (Math.abs(d) < tiny) d = tiny;
      c = 1 + aa / c; if (Math.abs(c) < tiny) c = tiny;
      d = 1 / d; const delta = d * c; h *= delta;
      if (Math.abs(delta - 1) < 2e-14) break;
    }
    return h;
  }
  function regularizedBeta(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    const bt = Math.exp(logGamma(a + b) - logGamma(a) - logGamma(b) + a * Math.log(x) + b * Math.log1p(-x));
    return clamp01(x < (a + 1) / (a + b + 2) ? bt * betaFraction(a, b, x) / a : 1 - bt * betaFraction(b, a, 1 - x) / b);
  }
  function gammaQ(a, x) {
    if (x <= 0) return 1;
    if (!Number.isFinite(x)) return 0;
    const scale = Math.exp(a * Math.log(x) - x - logGamma(a));
    if (x < a + 1) {
      let ap = a, sum = 1 / a, term = sum;
      for (let i = 1; i <= 1000; i++) { ap++; term *= x / ap; sum += term; if (Math.abs(term) < Math.abs(sum) * 2e-15) break; }
      return clamp01(1 - sum * scale);
    }
    let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
    for (let i = 1; i <= 1000; i++) {
      const an = -i * (i - a); b += 2;
      d = an * d + b; if (Math.abs(d) < 1e-300) d = 1e-300;
      c = b + an / c; if (Math.abs(c) < 1e-300) c = 1e-300;
      d = 1 / d; const delta = d * c; h *= delta;
      if (Math.abs(delta - 1) < 2e-15) break;
    }
    return clamp01(scale * h);
  }
  const normalPDF = x => Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI);
  const normalTwoSidedP = z => gammaQ(.5, z * z / 2);
  const tPDF = (t, df) => Math.exp(logGamma((df + 1) / 2) - logGamma(df / 2) - .5 * Math.log(df * Math.PI) - (df + 1) / 2 * Math.log1p(t * t / df));
  const tTwoSidedP = (t, df) => {
    const square = t * t;
    // Use the small argument directly near zero, avoiding 1 - x cancellation.
    return square < df ? 1 - regularizedBeta(square / (df + square), .5, df / 2) : regularizedBeta(df / (df + square), df / 2, .5);
  };
  const chiPDF = (x, df) => x <= 0 ? 0 : Math.exp((df / 2 - 1) * Math.log(x) - x / 2 - df / 2 * Math.log(2) - logGamma(df / 2));
  const chiP = (x, df) => gammaQ(df / 2, x / 2);
  function inverseTail(p, sf) {
    let lo = 0, hi = 1;
    while (sf(hi) > p && hi < 1e10) hi *= 2;
    for (let i = 0; i < 85; i++) { const mid = (lo + hi) / 2; if (sf(mid) > p) lo = mid; else hi = mid; }
    return (lo + hi) / 2;
  }
  function meanModel(id, p) {
    const df = p.n - 1, se = p.sd / Math.sqrt(p.n), statistic = (p.mean - p.mu) / se, isT = id === 't-statistic';
    const tail = isT ? x => tTwoSidedP(x, df) : normalTwoSidedP;
    const pv = tail(Math.abs(statistic)), critical = inverseTail(p.alpha, tail);
    return {se, statistic, df: isT ? df : null, p: pv, critical, lower: p.mean - critical * se, upper: p.mean + critical * se, reject: pv < p.alpha};
  }
  function chiModel(counts, alpha = .05) {
    const n = counts.reduce((s, x) => s + x, 0), expected = n / counts.length, df = counts.length - 1;
    const contributions = counts.map(x => expected ? (x - expected) ** 2 / expected : 0), statistic = contributions.reduce((s, x) => s + x, 0);
    const p = n ? chiP(statistic, df) : NaN, critical = inverseTail(alpha, x => chiP(x, df));
    return {n, expected, df, contributions, statistic, p, critical, adequate: expected >= 5, reject: p < alpha};
  }

  function meanCard(isT) {
    const symbol = isT ? 't' : 'Z', id = isT ? 't-statistic' : 'z-statistic';
    return {
      id, kind: 'application', year: null, yearLabel: 'Практика', author: 'Современный учебный разбор', group: 'inference',
      title: isT ? 't-статистика: когда разброс приходится оценивать' : 'Z-статистика: разница в единицах стандартной ошибки',
      shortTitle: isT ? 't-статистика' : 'Z-статистика',
      hook: isT ? 'Числа те же, а хвосты тяжелее. Почему меняется вывод?' : 'Среднее сдвинулось на два грамма. Это много или мало?',
      summary: isT ? 'Если разброс неизвестен, его оценивают по той же выборке. t-статистика учитывает эту дополнительную неопределённость через распределение с n−1 степенями свободы.' : 'Z-статистика показывает, на сколько стандартных ошибок выборочное среднее отстоит от проверяемого значения. При известном разбросе и нормальной модели её сравнивают со стандартным нормальным распределением.',
      problem: 'Мастерская фасует пакеты с пигментом по 100 г. В независимой выборке средний вес составил 102 г. Нужно проверить двустороннюю гипотезу H₀: μ = 100 г: нас интересует отклонение в любую сторону.',
      approach: isT ? 'Разброс s оцениваем по выбранным пакетам. Делим разницу среднего и нормы на s/√n; затем читаем хвосты t-распределения с n−1 степенями свободы.' : 'По отдельной надёжной калибровке известно σ = 4 г. Разброс среднего 16 независимых пакетов — σ/√16 = 1 г. Поэтому разнице в 2 г соответствует Z = 2.',
      prerequisites: isT ? ['student', 'z-statistic', 'sample-mean'] : ['gauss', 'sample-mean', 'neyman-pearson'],
      steps: [
        {title: 'Сначала задай вопрос', text: 'Нулевая гипотеза относится к среднему весу всей продукции: μ = 100 г. Альтернатива двусторонняя: μ ≠ 100 г. Уровень α выбирают до просмотра результата.'},
        {title: isT ? 'Оцени разброс по выборке' : 'Найди разброс среднего', text: isT ? 's² = Σ(Xᵢ−X̄)²/(n−1). Оценка s меняется от выборки к выборке. Стандартная ошибка среднего оценивается как s/√n.' : 'Стандартное отклонение σ описывает отдельные пакеты, а σ/√n — колебания среднего. Для σ = 4 г и n = 16 стандартная ошибка равна 1 г.'},
        {title: 'Раздели разницу на её стандартную ошибку', text: isT ? 'При X̄ = 102, μ₀ = 100, s = 4 и n = 16 получаем t = (102−100)/(4/√16) = 2. Число 2 не изменилось, но теперь закон сравнения — t₁₅.' : 'Z = (102−100)/1 = 2. Это две стандартные ошибки от нулевого среднего; знак показывает направление сдвига.'},
        {title: 'Посчитай площадь обоих хвостов', text: isT ? 'Для |t| ≥ 2 при df = 15 двустороннее p ≈ 0,06395. При α = 0,05 H₀ не отвергается. Порог около ±2,131: неопределённость оценки s делает его дальше от нуля.' : 'При H₀ стандартная нормальная модель даёт P(|Z| ≥ 2) ≈ 0,04550. При α = 0,05 результат попадает за порог ±1,960, и H₀ отвергается.'},
        {title: 'Проверь, что позволяет сказать модель', text: 'p-значение — вероятность не менее экстремальной статистики при H₀ и предпосылках модели. Оно не сообщает вероятность истинности H₀ и само по себе не измеряет важность разницы.'}
      ],
      formula: isT ? 'T = (X̄−μ₀)/(s/√n); при H₀ и независимой нормальной выборке T ~ tₙ₋₁.' : 'Z = (X̄−μ₀)/(σ/√n); при H₀ и независимой нормальной выборке с известным σ: Z ~ N(0,1).',
      example: isT ? 'Те же X̄ = 102 г, s = 4 г, n = 16: t = 2, df = 15, p ≈ 0,06395. 95%-интервал: [99,869; 104,131] г. Он включает 100 г. Замена известного σ выборочным s меняет распределение, даже когда знаменатели численно равны.' : 'При X̄ = 102 г, σ = 4 г и n = 16 получаем Z = 2, p ≈ 0,04550. 95%-интервал: [100,040; 103,960] г. Это вывод об общем среднем при модели, а не утверждение, что каждый пакет тяжелее нормы.',
      limitations: isT ? 'Точный закон tₙ₋₁ требует независимых нормальных наблюдений. При большом n метод часто используют приближённо; сильная асимметрия, выбросы и зависимость требуют отдельного внимания. Здесь одновыборочный тест: для двух независимых групп нужны другие стандартная ошибка и степени свободы.' : 'σ должна быть известна из обоснованной внешней модели, а не подставлена из этой же выборки. Для независимых нормальных наблюдений закон точен; для других распределений нормальное приближение требует условий. Большой объём не исправляет зависимость или систематический сдвиг.',
      takeaway: isT ? 'Если разброс оценён по данным, нужно учитывать случайность знаменателя. Одного значения t недостаточно: нужны ещё степени свободы.' : 'Сначала уточни, чей разброс стоит в знаменателе: одного наблюдения или среднего. Проверка среднего использует стандартную ошибку σ/√n.',
      applications: [{title: 'Контроль среднего веса', text: 'Сравниваем среднее независимых измерений с заранее заданной нормой. Чтобы понять практическую важность, вместе с p смотрим разницу в граммах и интервал.'}],
      sources: isT ? [nistMean, nistT] : [nistNormal, nistMean],
      quiz: isT ? {question: 'Z = 2 и t = 2 при df = 4. У какой статистики больше двустороннее p-значение?', options: ['У t: его хвосты тяжелее', 'У Z: нормальная кривая выше в центре', 'Они равны, потому что значения статистик совпали'], correct: 0, explanation: 'p зависит от закона сравнения. При df = 4 у t более тяжёлые хвосты: t = 2 даёт p около 0,116, а Z = 2 — около 0,0455.'} : {question: 'σ = 4 г. Что стоит в знаменателе Z при проверке среднего 16 независимых измерений?', options: ['4 г: разброс отдельного измерения', '1 г: стандартная ошибка среднего', '16: размер выборки'], correct: 1, explanation: 'Стандартная ошибка среднего равна σ/√n = 4/4 = 1 г. Именно с ней сравнивают отклонение среднего от μ₀.'},
      python: {filename: isT ? 't_statistic.py' : 'z_statistic.py', code: ''}
    };
  }

  function installData() {
    if (!global.ATLAS_DATA) return;
    const data = global.ATLAS_DATA;
    for (const isT of [false, true]) {
      const card = meanCard(isT);
      if (!data.ideas.some(x => x.id === card.id)) data.ideas.push(card);
      global.ATLAS_DIGEST ||= {};
      global.ATLAS_DIGEST[card.id] = {question: card.hook, idea: card.summary, example: card.example, carry: card.takeaway, context: card.problem};
      global.ATLAS_STORIES ||= {};
      global.ATLAS_STORIES[card.id] = {opening: card.hook};
    }
    const relations = [
      {from: 'sample-mean', to: 'z-statistic', type: 'prerequisite', label: 'От разброса среднего к статистике', explanation: 'Стандартная ошибка задаёт единицу, в которой измеряется отклонение среднего.'},
      {from: 'gauss', to: 'z-statistic', type: 'application', label: 'Нормальная модель в проверке среднего', explanation: 'При известном σ нормальное распределение задаёт хвостовые вероятности Z.'},
      {from: 'student', to: 't-statistic', type: 'application', label: 'От закона к наблюдаемой статистике', explanation: 'Оцениваем σ по той же выборке и сравниваем отношение с t-распределением.'},
      {from: 'z-statistic', to: 't-statistic', type: 'generalizes', label: 'Что меняет оценивание разброса', explanation: 'При одинаковой численной разнице t может дать большее p: его закон учитывает случайность s.'},
      {from: 't-statistic', to: 'neyman', type: 'application', label: 'Тест и интервал для одного среднего', explanation: 'Двусторонний t-тест и соответствующий интервал выражают одно правило через порог и через диапазон средних.'}
    ];
    for (const relation of relations) if (!data.relations.some(x => x.from === relation.from && x.to === relation.to && x.type === relation.type)) data.relations.push(relation);
    const pearson = data.ideas.find(x => x.id === 'pearson');
    if (pearson && !pearson.sources.some(s => s.url === nistChi.url)) pearson.sources.push(nistChi);
    global.ATLAS_FORMULAS ||= {};
    global.ATLAS_FORMULAS['z-statistic'] = {markdown: raw`$$Z=\frac{\overline X-\mu_0}{\sigma/\sqrt n},\qquad Z\mid H_0\sim\mathcal N(0,1)$$
$$p=2\bigl[1-\Phi(|z_{\mathrm{obs}}|)\bigr]$$`, explanation: 'σ известно; наблюдения независимы и нормальны. p относится к двусторонней альтернативе. Статистика — случайная величина, z_obs — её значение по полученной выборке.'};
    global.ATLAS_FORMULAS['t-statistic'] = {markdown: raw`$$T=\frac{\overline X-\mu_0}{s/\sqrt n},\qquad T\mid H_0\sim t_{n-1}$$
$$s^2=\frac{\sum_{i=1}^{n}(X_i-\overline X)^2}{n-1},\qquad p=2\bigl[1-F_{t_{n-1}}(|t_{\mathrm{obs}}|)\bigr]$$`, explanation: 's оценено по той же выборке. Закон точен при независимых нормальных наблюдениях; число степеней свободы равно n−1. Для двух выборок формула меняется.'};
    global.ATLAS_FORMULAS.pearson = {markdown: raw`$$\chi^2_{\mathrm{obs}}=\sum_{i=1}^{m}\frac{(O_i-E_i)^2}{E_i},\qquad \nu=m-1$$
$$p=P\!\left(\chi^2_{\nu}\ge\chi^2_{\mathrm{obs}}\right)$$`, explanation: 'Здесь m категорий с заранее заданными вероятностями, без оценивания параметров по этой выборке. Для честной кости m = 6, E_i = N/6, ν = 5. Закон сравнения приближённый; маленькие ожидаемые частоты ухудшают его точность.'};
  }

  const defaults = id => id === 'pearson' ? {counts: [9, 15, 10, 16, 7, 3], alpha: .05} : {mean: 102, mu: 100, sd: 4, n: 16, alpha: .05};
  function parameters(id) {
    const p = retained.get(id) || defaults(id);
    return {...p, ...(p.counts ? {counts: [...p.counts]} : {})};
  }
  function field(name, label, value, min, max, step, note = '') {
    return `<label class="stat-field"><span>${label}</span><input name="${name}" type="number" value="${value}" min="${min}" max="${max}" step="${step}" inputmode="decimal">${note ? `<small>${note}</small>` : ''}</label>`;
  }
  function alphaField(alpha) {
    return `<label class="stat-field"><span>Уровень α</span><select name="alpha">${[.10, .05, .01].map(x => `<option value="${x}"${alpha === x ? ' selected' : ''}>${format(x, 2)}</option>`).join('')}</select><small>Порог выбирают до анализа данных.</small></label>`;
  }
  function results(items) {
    return `<dl class="stat-results">${items.map(([label, value]) => `<div class="stat-result"><dt>${label}</dt><dd>${esc(value)}</dd></div>`).join('')}</dl>`;
  }
  function densityPlot(pdf, m, isT) {
    const limit = Math.max(4.5, Math.min(9, Math.max(Math.abs(m.statistic) * 1.15, m.critical * 1.2)));
    const left = 40, right = 610, top = 20, bottom = 246, ymax = .44;
    const X = x => left + (x + limit) / (2 * limit) * (right - left), Y = y => bottom - y / ymax * (bottom - top);
    const curve = fn => Array.from({length: 361}, (_, j) => { const x = -limit + 2 * limit * j / 360; return `${j ? 'L' : 'M'}${X(x).toFixed(2)},${Y(fn(x)).toFixed(2)}`; }).join(' ');
    const fillTail = (lo, hi) => {
      if (hi <= lo) return '';
      const pts = Array.from({length: 121}, (_, j) => { const x = lo + (hi - lo) * j / 120; return `L${X(x).toFixed(2)},${Y(pdf(x)).toFixed(2)}`; }).join(' ');
      return `<path d="M${X(lo)},${bottom} ${pts} L${X(hi)},${bottom} Z" fill="var(--reader-chart-secondary,#bc8754)" opacity=".25"/>`;
    };
    const a = Math.min(limit, Math.abs(m.statistic));
    let svg = `<svg class="stat-plot" viewBox="0 0 640 290" role="img" aria-label="${isT ? 't-распределение при нулевой гипотезе' : 'Стандартное нормальное распределение при нулевой гипотезе'}. Закрашены оба хвоста за модулем наблюдаемой статистики. p = ${esc(probability(m.p))}.">`;
    for (const x of [-Math.floor(limit), -2, 0, 2, Math.floor(limit)].filter((x,i,a) => a.indexOf(x) === i)) svg += `<line x1="${X(x)}" y1="${top}" x2="${X(x)}" y2="${bottom}" stroke="var(--reader-chart-grid,#ddd)"/><text x="${X(x)}" y="${bottom + 20}" text-anchor="middle">${x}</text>`;
    svg += fillTail(-limit, -a) + fillTail(a, limit);
    svg += `<path d="${curve(pdf)}" fill="none" stroke="var(--reader-chart-primary,#3470a6)" stroke-width="2.5"/>`;
    if (isT) svg += `<path d="${curve(normalPDF)}" fill="none" stroke="var(--reader-chart-muted,#67717d)" stroke-width="1.5" stroke-dasharray="5 6"/>`;
    for (const x of [-m.critical, m.critical]) if (Math.abs(x) < limit) svg += `<line x1="${X(x)}" y1="${top+16}" x2="${X(x)}" y2="${bottom}" stroke="var(--reader-chart-secondary,#bc8754)" stroke-dasharray="4 5"/>`;
    const statX = X(Math.max(-limit, Math.min(limit, m.statistic)));
    svg += `<line x1="${statX}" y1="${top}" x2="${statX}" y2="${bottom}" stroke="var(--ink,#223245)" stroke-width="1.5"/><text x="${Math.max(left + 20, Math.min(right - 20, statX))}" y="${top-5}" text-anchor="middle">${isT ? 't' : 'Z'} = ${format(m.statistic,2)}</text>`;
    svg += `<line x1="${left}" y1="${bottom}" x2="${right}" y2="${bottom}" stroke="var(--reader-chart-muted,#67717d)"/><text x="325" y="286" text-anchor="middle">Значение статистики при H₀</text></svg>`;
    if (m.critical > limit) svg += `<p class="visual-caption">Критические значения ±${format(m.critical,3)} находятся за пределами показанной шкалы.</p>`;
    return svg;
  }
  function meanHTML(id, p) {
    const isT = id === 't-statistic', symbol = isT ? 't' : 'Z', m = meanModel(id, p);
    const cross = isT ? normalTwoSidedP(Math.abs(m.statistic)) : tTwoSidedP(Math.abs(m.statistic), p.n - 1);
    return `${results([['Стандартная ошибка, г',format(m.se)], [symbol+' по выборке',format(m.statistic)], ['Двустороннее p',probability(m.p)], [isT?'Степени свободы':'Закон при H₀',isT?String(m.df):'N(0, 1)'], ['Критические значения','±'+format(m.critical,3)], [`${format((1-p.alpha)*100,0)}%-интервал μ, г`,`${format(m.lower,3)} … ${format(m.upper,3)}`]])}
      ${densityPlot(x => isT ? tPDF(x, m.df) : normalPDF(x), m, isT)}
      <div class="stat-key"><span>${isT ? 't, ν = '+m.df : 'N(0, 1)'}</span>${isT?'<span class="stat-normal">N(0, 1) для сравнения</span>':''}<span class="stat-tail">Площадь = p</span><span class="stat-critical">Пороги для α</span></div>
      <p class="stat-explanation">При α = ${format(p.alpha,2)} ${m.reject ? 'H₀ отвергается' : 'оснований отвергнуть H₀ недостаточно'}: p ${m.reject ? '<' : '≥'} α. ${m.reject ? 'Разницу среднего с нормой трудно объяснить одним выборочным разбросом в этой модели.' : 'Это не доказывает, что средний вес равен норме.'} Разница в исходных единицах — ${format(p.mean-p.mu,3)} г.</p>
      <p class="visual-caption">Закрашенные хвосты начинаются от ±|${symbol}|. Пунктирные пороги задаются α; это две разные границы. ${isT ? 'Если бы разброс был известен и мы использовали Z с тем же числом '+format(m.statistic,2)+', p было бы '+probability(cross)+'.' : 'Если те же '+format(p.sd)+' г оценены по этой выборке, нужен t с '+(p.n-1)+' степенями свободы: p = '+probability(cross)+'.'}</p>
      <p class="stat-warning visual-caption">${isT ? 'Точный t-закон требует независимой нормальной выборки. Выбор t сам по себе не устраняет выбросы или зависимость.' : 'Здесь σ считается известным по внешней калибровке. Подстановка s из этой же выборки меняет закон сравнения.'}</p>`;
  }
  function chiBars(p, m) {
    const left = 40, right = 610, top = 25, bottom = 195, max = Math.max(...p.counts,m.expected,1)*1.18;
    const X = i => left + (i+.5)/6*(right-left), Y = x => bottom - x/max*(bottom-top), barWidth = 38;
    let out = `<svg class="stat-plot" viewBox="0 0 640 245" role="img" aria-label="Число выпадений каждой из шести граней и одинаковое ожидаемое число ${format(m.expected)} при честной кости.">`;
    for(let i=0;i<6;i++) out += `<rect x="${X(i)-barWidth/2}" y="${Y(p.counts[i])}" width="${barWidth}" height="${bottom-Y(p.counts[i])}" fill="var(--reader-chart-primary,#3470a6)" opacity=".7"/><text x="${X(i)}" y="${Y(p.counts[i])-8}" text-anchor="middle">${p.counts[i]}</text><text x="${X(i)}" y="${bottom+21}" text-anchor="middle">${i+1}</text>`;
    out += `<line x1="${left}" y1="${Y(m.expected)}" x2="${right}" y2="${Y(m.expected)}" stroke="var(--reader-chart-secondary,#bc8754)" stroke-width="2" stroke-dasharray="5 4"/><line x1="${left}" y1="${bottom}" x2="${right}" y2="${bottom}" stroke="var(--reader-chart-muted,#67717d)"/><text x="325" y="241" text-anchor="middle">Грань кости · пунктир: ожидается N/6</text></svg>`;
    return out;
  }
  function chiTail(m) {
    const xmax = Math.max(20, Math.min(50,m.statistic*1.15)), X = x => 40+x/xmax*570, Y = y => 190-y/.17*170;
    const curve = Array.from({length:301},(_,j)=>{const x=xmax*j/300;return `${j?'L':'M'}${X(x)},${Y(chiPDF(x,m.df))}`;}).join(' ');
    const start = Math.min(m.statistic,xmax), tail=Array.from({length:151},(_,j)=>{const x=start+(xmax-start)*j/150;return `L${X(x)},${Y(chiPDF(x,m.df))}`;}).join(' ');
    return `<svg class="stat-plot" viewBox="0 0 640 238" role="img" aria-label="Распределение хи-квадрат с пятью степенями свободы. Верхний хвост за наблюдаемой статистикой имеет площадь p = ${esc(probability(m.p))}."><path d="M${X(start)},190 ${tail} L610,190 Z" fill="var(--reader-chart-secondary,#bc8754)" opacity=".25"/><path d="${curve}" fill="none" stroke="var(--reader-chart-primary,#3470a6)" stroke-width="2.5"/><line x1="${X(m.critical)}" y1="20" x2="${X(m.critical)}" y2="190" stroke="var(--reader-chart-secondary,#bc8754)" stroke-dasharray="4 5"/><line x1="${X(start)}" y1="20" x2="${X(start)}" y2="190" stroke="var(--ink,#223245)"/><line x1="40" y1="190" x2="610" y2="190" stroke="var(--reader-chart-muted,#67717d)"/>${[0,5,10,15,20,...(xmax>25?[30,40]:[])].filter(x=>x<=xmax).map(x=>`<text x="${X(x)}" y="210" text-anchor="middle">${x}</text>`).join('')}<text x="325" y="235" text-anchor="middle">χ² при H₀ · ν = 5</text></svg>`;
  }
  function chiHTML(p) {
    const m = chiModel(p.counts,p.alpha);
    if (!m.n) return '<p class="stat-error" role="alert">Добавь хотя бы одно выпадение: при N = 0 ожидаемые частоты и критерий не определены.</p>';
    return `${results([['Всего бросков N',String(m.n)],['Ожидается на грань',format(m.expected,3)],['Наблюдаемое χ²',format(m.statistic)],['Степени свободы',String(m.df)],['Приближённое p',probability(m.p)],['Порог для α',format(m.critical,3)]])}${chiBars(p,m)}
      <div class="stat-table-wrap"><table><caption>Из чего складывается χ²</caption><thead><tr><th scope="col">Грань</th><th scope="col">Oᵢ</th><th scope="col">Eᵢ = N/6</th><th scope="col">(Oᵢ−Eᵢ)²/Eᵢ</th></tr></thead><tbody>${p.counts.map((x,i)=>`<tr><th scope="row">${i+1}</th><td>${x}</td><td>${format(m.expected,3)}</td><td>${format(m.contributions[i],3)}</td></tr>`).join('')}</tbody></table></div>
      ${chiTail(m)}<div class="stat-key"><span>χ², ν = 5</span><span class="stat-tail">Верхний хвост = p</span><span class="stat-critical">Порог для α</span></div>
      ${m.adequate ? `<p class="stat-explanation">При α = ${format(p.alpha,2)} ${m.reject ? 'модель равновероятных граней отвергается' : 'оснований отвергнуть модель равновероятных граней недостаточно'}. p = ${probability(m.p)} относится к независимым броскам и вероятностям 1/6. ${m.reject ? 'Причину расхождения этот расчёт не устанавливает.' : 'Такой результат не доказывает, что кость честная.'}</p>` : `<p class="stat-warning">Ожидаемая частота ${format(m.expected,3)} меньше 5. Показанное p — приближение, которому здесь не стоит доверять для решения. Для малого N нужен точный мультиномиальный расчёт или калиброванный случайный эксперимент.</p>`}
      <p class="visual-caption">Шесть частот связаны условием, что их сумма равна N, поэтому ν = 6−1 = 5. Вероятности 1/6 заданы заранее; параметры по этим броскам не оценивались. Для таблицы сопряжённости степени свободы считаются иначе.${m.statistic>50?' Наблюдаемое χ² находится правее показанной шкалы.':''}</p>`;
  }
  function mount(container,id) {
    unmount(); if (!container) return;
    const p = parameters(id), isChi=id==='pearson', isT=id==='t-statistic';
    const title=isChi?'Шесть граней: от частот к хвостовой вероятности':isT?'Та же разница, неизвестный разброс':'В каких единицах измеряется отклонение среднего?';
    const controls=isChi ? p.counts.map((x,i)=>field('count-'+i,'Выпадения грани '+(i+1),x,0,10000,1)).join('') : field('mean','Выборочное среднее X̄, г',p.mean,-10000,10000,'any')+field('mu','Проверяемое среднее μ₀, г',p.mu,-10000,10000,'any')+field('sd',isT?'Выборочный разброс s, г':'Известный разброс σ, г',p.sd,.01,10000,'any',isT?'Оценка по этой же выборке.':'Задан по внешней калибровке.')+field('n','Независимых измерений n',p.n,2,10000,1);
    container.innerHTML=`<section class="stat-lab"><h3>${title}</h3><p class="visual-caption">${isChi?'Нулевая модель: у каждой грани вероятность 1/6. Измени целые количества выпадений; общее число бросков пересчитается.':'Нулевая гипотеза: μ = μ₀. Альтернатива двусторонняя. Меняй по одному числу и следи за стандартной ошибкой, статистикой и p.'}</p><div class="stat-controls${isChi?' stat-die-controls':''}">${controls}${alphaField(p.alpha)}</div><div class="stat-actions">${isChi?'<button type="button" data-stat-preset="equal">Поровну: χ² = 0</button><button type="button" data-stat-preset="unequal">60 бросков с расхождением</button>':'<button type="button" data-stat-preset="default">Вернуть пример: 102 г</button><button type="button" data-stat-preset="null">Среднее равно норме</button>'}</div><div class="stat-output" aria-live="polite" aria-atomic="false"></div></section>`;
    const output=container.querySelector('.stat-output');
    const read=()=>{
      const names=isChi?[...p.counts.map((_,i)=>'count-'+i),'alpha']:['mean','mu','sd','n','alpha'];
      const values={};
      for(const name of names){const input=container.querySelector(`[name="${name}"]`),value=Number(input.value);if(input.value===''||!Number.isFinite(value)||!input.checkValidity()){output.innerHTML='<p class="stat-error" role="alert">Проверь числа: n — целое от 2 до 10 000, разброс положительный, частоты — целые от 0 до 10 000.</p>';return;}values[name]=value;}
      if(isChi){p.counts=p.counts.map((_,i)=>values['count-'+i]);p.alpha=values.alpha;}else Object.assign(p,values);
      retained.set(id,parametersCopy(p));output.innerHTML=isChi?chiHTML(p):meanHTML(id,p);
    };
    const inputHandler=()=>read();
    const clickHandler=e=>{const button=e.target.closest('[data-stat-preset]');if(!button)return;const preset=button.dataset.statPreset;
      if(isChi){const values=preset==='equal'?[10,10,10,10,10,10]:[9,15,10,16,7,3];values.forEach((v,i)=>container.querySelector(`[name="count-${i}"]`).value=v);}
      else if(preset==='null')container.querySelector('[name="mean"]').value=container.querySelector('[name="mu"]').value;
      else{const d=defaults(id);for(const name of ['mean','mu','sd','n','alpha'])container.querySelector(`[name="${name}"]`).value=d[name];}read();};
    container.addEventListener('input',inputHandler);container.addEventListener('change',inputHandler);container.addEventListener('click',clickHandler);
    current={container,inputHandler,clickHandler};read();
  }
  const parametersCopy=p=>({...p,...(p.counts?{counts:[...p.counts]}:{})});
  function unmount(){if(current){const c=current;c.container.removeEventListener('input',c.inputHandler);c.container.removeEventListener('change',c.inputHandler);c.container.removeEventListener('click',c.clickHandler);current=null;}}

  function art(id,variant='paper') {
    const isT=id==='t-statistic',pdf=isT?x=>tPDF(x,15):normalPDF, limit=4, X=x=>95+(x+limit)/(2*limit)*710,Y=y=>430-y*610;
    const curve=fn=>Array.from({length:201},(_,j)=>{const x=-limit+2*limit*j/200;return `${j?'L':'M'}${X(x)},${Y(fn(x))}`;}).join(' ');
    const a=2,fill=(lo,hi)=>`M${X(lo)},430 `+Array.from({length:101},(_,j)=>{const x=lo+(hi-lo)*j/100;return `L${X(x)},${Y(pdf(x))}`;}).join(' ')+` L${X(hi)},430 Z`;
    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 900 600" role="img" aria-label="${isT?'t-статистика: оценённый разброс, более тяжёлые хвосты':'Z-статистика: отклонение среднего в стандартных ошибках'}" class="atlas-art atlas-art--${esc(variant)}" style="display:block;width:100%;height:auto;overflow:visible"><text x="75" y="82" fill="var(--ink,#263242)" font-size="24">${isT?'РАЗБРОС ОЦЕНЁН ПО ВЫБОРКЕ':'РАЗНИЦА / СТАНДАРТНАЯ ОШИБКА'}</text><text x="75" y="120" fill="var(--muted,#67717d)" font-size="20">${isT?'t = 2 · ν = 15 · p ≈ 0,06395':'Z = 2 · известное σ · p ≈ 0,0455'}</text><path d="${fill(-4,-a)} ${fill(a,4)}" fill="var(--reader-chart-secondary,#bc8754)" opacity=".3"/>${isT?`<path d="${curve(normalPDF)}" fill="none" stroke="var(--reader-chart-muted,#67717d)" stroke-width="2" stroke-dasharray="7 8"/>`:''}<path d="${curve(pdf)}" fill="none" stroke="var(--reader-chart-primary,#3470a6)" stroke-width="3"/><path d="M95,430 H805" stroke="var(--ink,#263242)" fill="none"/>${[-2,0,2].map(x=>`<path d="M${X(x)},430 v9" stroke="var(--ink,#263242)"/><text x="${X(x)}" y="475" text-anchor="middle" fill="var(--ink,#263242)" font-size="24">${x}</text>`).join('')}<text x="450" y="536" text-anchor="middle" fill="var(--muted,#67717d)" font-size="20">${isT?'Неопределённость s меняет закон сравнения':'Оба хвоста за ±2 дают двустороннее p'}</text></svg>`;
  }

  // A complete standard-library Python counterpart, populated with current controls.
  function pythonCode(id) {
    const p=parameters(id);
    const common=`import math\n\ndef beta_fraction(a, b, x):\n    tiny = 1e-300\n    c, d = 1.0, 1 - (a+b)*x/(a+1)\n    d = 1 / (d if abs(d) > tiny else tiny)\n    h = d\n    for m in range(1, 501):\n        for aa in (m*(b-m)*x/((a+2*m-1)*(a+2*m)), -(a+m)*(a+b+m)*x/((a+2*m)*(a+2*m+1))):\n            d = 1 + aa*d\n            c = 1 + aa/c\n            d = 1 / (d if abs(d) > tiny else tiny)\n            c = c if abs(c) > tiny else tiny\n            delta = d*c\n            h *= delta\n        if abs(delta-1) < 2e-14:\n            break\n    return h\n\ndef beta_i(x, a, b):\n    if x <= 0: return 0.0\n    if x >= 1: return 1.0\n    bt = math.exp(math.lgamma(a+b)-math.lgamma(a)-math.lgamma(b)+a*math.log(x)+b*math.log1p(-x))\n    return bt*beta_fraction(a,b,x)/a if x < (a+1)/(a+b+2) else 1-bt*beta_fraction(b,a,1-x)/b\n\ndef gamma_q(a, x):\n    if x <= 0: return 1.0\n    scale = math.exp(a*math.log(x)-x-math.lgamma(a))\n    if x < a+1:\n        ap, total, term = a, 1/a, 1/a\n        for _ in range(1000):\n            ap += 1\n            term *= x/ap\n            total += term\n            if abs(term) < abs(total)*2e-15: break\n        return max(0.0, min(1.0, 1-total*scale))\n    b, c = x+1-a, 1e300\n    d = 1/b\n    h = d\n    for i in range(1,1001):\n        an = -i*(i-a)\n        b += 2\n        d = an*d+b\n        c = b+an/c\n        d = 1/(d if abs(d)>1e-300 else 1e-300)\n        c = c if abs(c)>1e-300 else 1e-300\n        delta = d*c\n        h *= delta\n        if abs(delta-1)<2e-15: break\n    return max(0.0, min(1.0,scale*h))\n\ndef inverse_tail(p, sf):\n    lo, hi = 0.0, 1.0\n    while sf(hi)>p: hi *= 2\n    for _ in range(85):\n        mid = (lo+hi)/2\n        if sf(mid)>p: lo = mid\n        else: hi = mid\n    return (lo+hi)/2\n\n`;
    if(id==='pearson')return '# Python 3; только стандартная библиотека. Текущие параметры карточки.\n'+common+`counts = ${JSON.stringify(p.counts)}\nalpha = ${p.alpha}\nn = sum(counts)\nif not n: raise ValueError("Нужно хотя бы одно наблюдение")\nexpected = n/6\ndf = 5\nstatistic = sum((x-expected)**2/expected for x in counts)\np_value = gamma_q(df/2,statistic/2)\ncritical = inverse_tail(alpha,lambda x:gamma_q(df/2,x/2))\nprint(f"N={n}, E={expected:.6f}, chi2={statistic:.6f}, df={df}")\nprint(f"p={p_value:.9g}, critical={critical:.6f}")\nif expected < 5:\n    print("Малые ожидаемые частоты: приближение ненадёжно, решение не формулируем.")\nelse:\n    print("H0 отвергается" if p_value < alpha else "Недостаточно оснований отвергнуть H0")\n`;
    const isT=id==='t-statistic';
    return '# Python 3; только стандартная библиотека. Текущие параметры карточки.\n'+common+`mean, mu0, sd, n, alpha = ${p.mean}, ${p.mu}, ${p.sd}, ${p.n}, ${p.alpha}\n# ${isT?'sd — выборочное s; предполагается независимая нормальная выборка.':'sd — известное sigma из внешней модели; предполагается независимая нормальная выборка.'}\nse = sd/math.sqrt(n)\nstatistic = (mean-mu0)/se\ndf = n-1\nsf = ${isT?'lambda x: 1-beta_i(x*x/(df+x*x),.5,df/2) if x*x < df else beta_i(df/(df+x*x),df/2,.5)':'lambda x: math.erfc(x/math.sqrt(2))'}\np_value = sf(abs(statistic))\ncritical = inverse_tail(alpha,sf)\nprint(f"SE={se:.6f}, statistic={statistic:.6f}${isT?', df={df}':''}")\nprint(f"two-sided p={p_value:.9g}, critical=+/-{critical:.6f}")\nprint(f"{100*(1-alpha):g}%-interval: [{mean-critical*se:.6f}, {mean+critical*se:.6f}]")\nprint("H0 отвергается" if p_value < alpha else "Недостаточно оснований отвергнуть H0")\n`;
  }
  function codeFor(id){return {filename:id==='pearson'?'chi_square.py':id==='t-statistic'?'t_statistic.py':'z_statistic.py',code:pythonCode(id),command:'python3 '+(id==='pearson'?'chi_square.py':id==='t-statistic'?'t_statistic.py':'z_statistic.py'),kind:'statistics'};}

  global.AtlasStatistics={has:id=>handled.has(id),mount,unmount,getParams:parameters,codeFor,math:{logGamma,regularizedBeta,gammaQ,normalPDF,normalTwoSidedP,tPDF,tTwoSidedP,chiP,inverseTail,meanModel,chiModel}};
  installData();
  const base=global.IdeaVisuals;
  if(base){
    const originalMount=base.mount,originalUnmount=base.unmount,originalParams=base.getParams;
    base.configs['z-statistic']={kind:'atlas-z',title:'Z-статистика'};
    base.configs['t-statistic']={kind:'atlas-t',title:'t-статистика'};
    base.mount=(container,id)=>{if(handled.has(id)){originalUnmount();mount(container,id);}else{unmount();originalMount(container,id);}};
    base.unmount=()=>{unmount();originalUnmount();};
    base.getParams=id=>handled.has(id)?null:originalParams(id);
  }
  const baseArt=global.AtlasArt?.svg;
  if(baseArt)global.AtlasArt.svg=(id,variant)=>id==='z-statistic'||id==='t-statistic'?art(id,variant):baseArt(id,variant);
})(typeof window==='undefined'?globalThis:window);
