from DrissionPage import ChromiumPage, ChromiumOptions, errors
import time
import re

co = ChromiumOptions().set_paths(browser_path=r'/Applications/Google Chrome.app').set_local_port(9202)
# Windows example if needed:
# co = ChromiumOptions().set_paths(browser_path=r'C:\Program Files (x86)\Microsoft\Edge\Application\msedge.exe').set_local_port(9202)

page = ChromiumPage(co)

# Matches: "gnomAD v4.1.0 variants (1023)" (count may include commas)
VARIANT_RE = re.compile(r'gnomAD\s+v[\d.]+\s+variants\s*\(([\d,]+)\)', re.IGNORECASE)


def get_variant_count():
    """
    Fast + robust:
    find the single span (or any element) whose text contains:
      'gnomAD' and 'variants ('
    then regex-extract the number.
    """
    ele = None
    for _ in range(6):
        # Try span first (most likely)
        ele = page.ele('xpath://span[contains(., "gnomAD v4.1.0 variants (")]')
        if ele:
            break
        time.sleep(0.5)

    if not ele:
        return None

    text = (ele.text or "").strip()
    m = VARIANT_RE.search(text)
    if m:
        return int(m.group(1).replace(",", ""))

    # If version text differs slightly but still has 'variants (###)', fallback regex
    m2 = re.search(r'variants\s*\(([\d,]+)\)', text, re.IGNORECASE)
    if m2:
        return int(m2.group(1).replace(",", ""))

    return None


def crawl(gene_name):
    try:
        page.get('https://gnomad.broadinstitute.org/')
        time.sleep(1)  # down time

        #dropdown = page.ele('tag:select@@class=Select-sc-1lkyg9e-0 iPbHsk')
        #dropdown.select('gnomad_r4')

        page.ele('tag:input@@data-testid=searchbox-input').input(gene_name + '\n')
        time.sleep(100)  # down time

        variant_count_all = get_variant_count()

        page.ele('tag:input@@id=variant-consequence-category-filter-synonymous').click()
        page.ele('tag:input@@id=variant-consequence-category-filter-other').click()
        time.sleep(1)  # down time

        # Extract count AFTER filters (fast, element-based)
        variant_count = get_variant_count()

        page.ele('tag:button@@text()=Export variants to CSV').click()

        # Print gene + count
        print(f"{gene_name}\t{variant_count_all if variant_count_all is not None else 'NA'}\t{variant_count if variant_count is not None else 'NA'}")

    except errors.ElementNotFoundError:
        print(f"{gene_name}\tNA")


if __name__ == '__main__':
    with open('1.txt', 'r', encoding='utf-8') as file:
        content = file.read().replace('\n', '')
        string_list = content.split(',')

    for item in string_list:
        crawl(item.strip())

